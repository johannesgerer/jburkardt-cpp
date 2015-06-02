# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "compass_search.hpp"

int main ( );
double beale ( int m, double x[] );
void beale_test ( );
double bohach1 ( int m, double x[] );
void bohach1_test ( );
double bohach2 ( int m, double x[] );
void bohach2_test ( );
double broyden ( int m, double x[] );
void broyden_test ( );
double extended_rosenbrock ( int m, double x[] );
void extended_rosenbrock_test ( );
double goldstein_price ( int m, double x[] );
void goldstein_price_test ( );
double himmelblau ( int m, double x[] );
void himmelblau_test ( );
double local ( int m, double x[] );
void local_test ( );
double mckinnon ( int m, double x[] );
void mckinnon_test ( );
void mckinnon_parameters ( string action, double &phi, double &tau, double &theta );
double powell ( int m, double x[] );
void powell_test ( );
double rosenbrock ( int m, double x[] );
void rosenbrock_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for COMPASS_SEARCH_PRB.
//
//  Discussion:
//
//    COMPASS_SEARCH_PRB tests the COMPASS_SEARCH library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "COMPASS_SEARCH_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the COMPASS_SEARCH library.\n";

  beale_test ( );
  bohach1_test ( );
  bohach2_test ( );
  broyden_test ( );
  extended_rosenbrock_test ( );
  goldstein_price_test ( );
  himmelblau_test ( );
  local_test ( );
  mckinnon_test ( );
  powell_test ( );
  rosenbrock_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "COMPASS_SEARCH_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void beale_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BEALE_TEST works with the Beale function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "BEALE_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Beale function.\n";
  delta_tol = 0.00001;
  delta = 0.1;
  k_max = 20000;

  x0[0] = 1.0;
  x0[1] = 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << beale ( m, x0 ) << "\n";

  x = compass_search ( beale, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Repeat with more difficult start.
//
  x = new double[m];
  x0[0] = 1.0;
  x0[1] = 4.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << beale ( m, x0 ) << "\n";
  delete [] x;

  x = compass_search ( beale, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx 
       << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 3.0;
  x[1] = 0.5;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << beale ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void bohach1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BOHACH1_TEST works with the Bohachevsky function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "BOHACH1_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Bohachevsky function #1.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = 0.5;
  x0[1] = 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << bohach1 ( m, x0 ) << "\n";

  x = compass_search ( bohach1, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx
       << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 0.0;
  x[1] = 0.0;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << bohach1 ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void bohach2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BOHACH2_TEST works with the Bohachevsky function #2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "BOHACH2_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Bohachevsky function #2.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = 0.6;
  x0[1] = 1.3;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << bohach2 ( m, x0 ) << "\n";

  x = compass_search ( bohach2, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 0.0;
  x[1] = 0.0;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << bohach2 ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void broyden_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BROYDEN_TEST works with the Broyden function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "BROYDEN_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Broyden function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = - 0.9;
  x0[1] = - 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << broyden ( m, x0 ) << "\n";

  x = compass_search ( broyden, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = -0.511547141775014;
  x[1] = -0.398160951772036;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << broyden ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void extended_rosenbrock_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    EXTENDED_ROSENBROCK_TEST works with the extended Rosenbrock function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 4;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "EXTENDED_ROSENBROCK_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the extended Rosenbrock function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = -1.2;
  x0[1] =  1.0;
  x0[2] = -1.5;
  x0[3] =  1.2;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << extended_rosenbrock ( m, x0 ) << "\n";

  x = compass_search ( extended_rosenbrock, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 1.0;
  x[1] = 1.0;
  x[2] = 1.0;
  x[3] = 1.0;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << extended_rosenbrock ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void goldstein_price_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    GOLDSTEIN_PRICE_TEST works with the Goldstein-Price function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "GOLDSTEIN_PRICE_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Goldstein-Price function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = -0.5;
  x0[1] =  0.25;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << goldstein_price ( m, x0 ) << "\n";

  x = compass_search ( goldstein_price, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Repeat with more difficult start.
//
  x0[0] = -4.0;
  x0[1] = 5.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << goldstein_price ( m, x0 ) << "\n";

  x = compass_search ( goldstein_price, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 0.0;
  x[1] = -1.0;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << goldstein_price ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void himmelblau_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    HIMMELBLAU_TEST works with the Himmelblau function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "HIMMELBLAU_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Himmelblau function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = 1.0;
  x0[1] = 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << himmelblau ( m, x0 ) << "\n";

  x = compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x0[0] = -1.0;
  x0[1] = 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << himmelblau ( m, x0 ) << "\n";

  x = compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x0[0] = -1.0;
  x0[1] = -1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << himmelblau ( m, x0 ) << "\n";

  x = compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x0[0] = 1.0;
  x0[1] = -1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << himmelblau ( m, x0 ) << "\n";

  x = compass_search ( himmelblau, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate Himmelblau minimizers.
//
  x = new double[m];
  x[0] = 3.0;
  x[1] = 2.0;
  r8vec_print ( m, x, "  Correct Himmelblau minimizer X1*:" );
  cout << "\n";
  cout << "  F(X*) = " << himmelblau ( m, x ) << "\n";
  delete [] x;

  x = new double[m];
  x[0] = 3.58439;
  x[1] = -1.84813;
  r8vec_print ( m, x, "  Correct Himmelblau minimizer X2*:" );
  cout << "\n";
  cout << "  F(X*) = " << himmelblau ( m, x ) << "\n";
  delete [] x;

  x = new double[m];
  x[0] = -3.77934;
  x[1] = -3.28317;
  r8vec_print ( m, x, "  Correct Himmelblau minimizer X3*:" );
  cout << "\n";
  cout << "  F(X*) = " << himmelblau ( m, x ) << "\n";
  delete [] x;

  x = new double[m];
  x[0] = -2.80512;
  x[1] = 3.13134;
  r8vec_print ( m, x, "  Correct Himmelblau minimizer X4*:" );
  cout << "\n";
  cout << "  F(X*) = " << himmelblau ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void local_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_TEST works with the Local function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "LOCAL_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Local function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = 1.0;
  x0[1] = 1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << local ( m, x0 ) << "\n";

  x = compass_search ( local, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate local minimizer.
//
  x = new double[m];
  x[0] = 0.2858054412;
  x[1] = 0.2793263206;
  r8vec_print ( m, x, "  Correct local minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << local ( m, x ) << "\n";
  delete [] x;
//
//  Try for global minimizer.
//
  x0[0] = -15.0;
  x0[1] = -35.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << local ( m, x0 ) << "\n";

  x = compass_search ( local, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate global minimizer.
//
  x = new double[m];
  x[0] = -21.02667179;
  x[1] = -36.75997872;
  r8vec_print ( m, x, "  Correct global minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << local ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void mckinnon_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    MCKINNON_TEST works with the McKinnon function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double phi;
  double tau;
  double theta;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "MCKINNON_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the McKinnon function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;
//
//  Test 1
//
  a = ( 1.0 + sqrt ( 33.0 ) ) / 8.0;
  b = ( 1.0 - sqrt ( 33.0 ) ) / 8.0;

  phi = 10.0;
  tau = 1.0;
  theta = 15.0;

  mckinnon_parameters ( "set", phi, tau, theta );

  x0[0] = a;
  x0[1] = b;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "  PHI = " << phi << " TAU = " << tau << " THETA = " << theta << "\n";
  cout << "\n";
  cout << "  F(X0) = " << mckinnon ( m, x0 ) << "\n";

  x = compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x = new double[m];
  x[0] = 0.0;
  x[1] = -0.5;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << mckinnon ( m, x ) << "\n";
  delete [] x;
//
//  Test 2
//
  a = ( 1.0 + sqrt ( 33.0 ) ) / 8.0;
  b = ( 1.0 - sqrt ( 33.0 ) ) / 8.0;

  phi = 60.0;
  tau = 2.0;
  theta = 6.0;

  mckinnon_parameters ( "set", phi, tau, theta );

  x0[0] = a;
  x0[1] = b;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "  PHI = " << phi << " TAU = " << tau << " THETA = " << theta << "\n";
  cout << "\n";
  cout << "  F(X0) = " << mckinnon ( m, x0 ) << "\n";

  x = compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x = new double[m];
  x[0] = 0.0;
  x[1] = -0.5;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << mckinnon ( m, x ) << "\n";
  delete [] x;
//
//  Test 3
//
  a = ( 1.0 + sqrt ( 33.0 ) ) / 8.0;
  b = ( 1.0 - sqrt ( 33.0 ) ) / 8.0;

  phi = 4000.0;
  tau = 3.0;
  theta = 6.0;

  mckinnon_parameters ( "set", phi, tau, theta );

  x0[0] = a;
  x0[1] = b;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "  PHI = " << phi << " TAU = " << tau << " THETA = " << theta << "\n";
  cout << "\n";
  cout << "  F(X0) = " << mckinnon ( m, x0 ) << "\n";

  x = compass_search ( mckinnon, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;

  x = new double[m];
  x[0] = 0.0;
  x[1] = -0.5;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << mckinnon ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void powell_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    POWELL_TEST works with the Powell function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 4;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "POWELL_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Powell function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] =  3.0;
  x0[1] = -1.0;
  x0[2] =  0.0;
  x0[3] =  1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << powell ( m, x0 ) << "\n";

  x = compass_search ( powell, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << powell ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

void rosenbrock_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK_TEST works with the Rosenbrock function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double delta;
  double delta_tol;
  double fx;
  int k;
  int k_max;
  int m = 2;
  double *x;
  double *x0;

  x0 = new double[m];
  cout << "\n";
  cout << "ROSENBROCK_TEST:\n";
  cout << "  Test COMPASS_SEARCH with the Rosenbrock function.\n";
  delta_tol = 0.00001;
  delta = 0.3;
  k_max = 20000;

  x0[0] = -1.2;
  x0[1] =  1.0;
  r8vec_print ( m, x0, "  Initial point X0:" );
  cout << "\n";
  cout << "  F(X0) = " << rosenbrock ( m, x0 ) << "\n";

  x = compass_search ( rosenbrock, m, x0, delta_tol, delta, k_max, fx, k );
  r8vec_print ( m, x, "  Estimated minimizer X1:" );
  cout << "\n";
  cout << "  F(X1) = " << fx << " number of steps = " << k << "\n";
  delete [] x;
//
//  Demonstrate correct minimizer.
//
  x = new double[m];
  x[0] = 1.0;
  x[1] = 1.0;

  r8vec_print ( m, x, "  Correct minimizer X*:" );
  cout << "\n";
  cout << "  F(X*) = " << rosenbrock ( m, x ) << "\n";
  delete [] x;

  delete [] x0;

  return;
}
//****************************************************************************80

double beale ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    BEALE computes the Beale function.
//
//  Discussion:
//
//    This function has a global minimizer:
//
//      X = ( 3.0, 0.5 )
//
//    for which
//
//      F(X) = 0.
//
//    For a relatively easy computation, start the search at
//
//      X = ( 1.0, 1.0 )
//
//    A harder computation starts at
//
//      X = ( 1.0, 4.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Evelyn Beale,
//    On an Iterative Method for Finding a Local Minimum of a Function
//    of More than One Variable,
//    Technical Report 25, 
//    Statistical Techniques Research Group,
//    Princeton University, 1958.
//
//    Richard Brent,
//    Algorithms for Minimization with Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double BEALE, the value of the function at X.
//
{
  double f;
  double f1;
  double f2;
  double f3;

  f1 = 1.5   - x[0] * ( 1.0 - x[1]    );
  f2 = 2.25  - x[0] * ( 1.0 - pow ( x[1], 2 ) );
  f3 = 2.625 - x[0] * ( 1.0 - pow ( x[1], 3 ) );

  f = f1 * f1 + f2 * f2 + f3 * f3;
 
  return f;
}
//****************************************************************************80

double bohach1 ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOHACH1 evaluates the Bohachevsky function #1.
//
//  Discussion:
//
//    The minimizer is
//
//      X* = [ 0.0, 0.0 ]
//      F(X*) = 0.0
//
//    Suggested starting point:
//
//      X = [ 0.5, 1.0 ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer Verlag, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double BOHACH1, the value of the function at X.
//
{
  double f;
  double pi = 3.141592653589793;

  f =       x[0] * x[0] - 0.3 * cos ( 3.0 * pi * x[0] ) 
    + 2.0 * x[1] * x[1] - 0.4 * cos ( 4.0 * pi * x[1] ) 
    + 0.7;
 
  return f;
}
//****************************************************************************80

double bohach2 ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOHACH2 evaluates the Bohachevsky function #2.
//
//  Discussion:
//
//    The minimizer is
//
//      X* = [ 0.0, 0.0 ]
//      F(X*) = 0.0
//
//    Suggested starting point:
//
//      X = [ 0.6, 1.3 ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer Verlag, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double BOHACH2, the value of the function at X.
//
{
  double f;
  double pi = 3.141592653589793;

  f =       x[0] * x[0] 
    + 2.0 * x[1] * x[1] 
    - 0.3 * cos ( 3.0 * pi * x[0] ) * cos ( 4.0 * pi * x[1] ) 
    + 0.3;

  return f;
}
//****************************************************************************80

double broyden ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    BROYDEN computes the two dimensional modified Broyden function.
//
//  Discussion:
//
//    This function has a global minimizer:
//
//      X = ( -0.511547141775014, -0.398160951772036 )
//
//    for which
//
//      F(X) = 1.44E-04
//
//    Start the search at
//
//      X = ( -0.9, -1.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Broyden,
//    A class of methods for solving nonlinear simultaneous equations,
//    Mathematics of Computation,
//    Volume 19, 1965, pages 577-593.
//
//    Jorge More, Burton Garbow, Kenneth Hillstrom,
//    Testing unconstrained optimization software,
//    ACM Transactions on Mathematical Software,
//    Volume 7, Number 1, March 1981, pages 17-41. 
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double BROYDEN, the value of the function at X.
//
{
  double f;
  double f1;
  double f2;
  double p;

  f1 = r8_abs ( ( 3.0 -           x[0] ) * x[0] - 2.0 * x[1] + 1.0 );
  f2 = r8_abs ( ( 3.0 - 2.0 * x[1] ) * x[1] -           x[0] + 1.0 );

  p = 3.0 / 7.0;

  f = pow ( f1, p ) + pow ( f2, p );
 
  return f;
}
//****************************************************************************80

double extended_rosenbrock ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXTENDED_ROSENBROCK computes the extended Rosenbrock function.
//
//  Discussion:
//
//    The number of dimensions is arbitrary, except that it must be even.
//
//    There is a global minimum at X* = (1,1,&), F(X*) = 0.
//
//    The contours are sharply twisted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Rosenbrock,
//    An Automatic Method for Finding the Greatest or Least Value of a Function,
//    Computer Journal,
//    Volume 3, 1960, pages 175-184.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double EXTENDED_ROSENBROCK, the value of the function at X.
//
{
  double f;
  double f1;
  double f2;
  int i;

  f = 0.0;

  for ( i = 0; i < m - 1; i = i + 2 )
  {
    f1 = 1.0 - x[i];
    f2 = 10.0 * ( x[i+1] - x[i] * x[i] );
    f = f + f1 * f1 + f2 * f2;
  }

  return f;
}
//****************************************************************************80

double goldstein_price ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    GOLDSTEIN_PRICE evaluates the Goldstein-Price polynomial.
//
//  Discussion:
//
//    The minimizer is
//
//      X* = [ 0.0, -1.0 ]
//      F(X*) = 3.0
//
//    Suggested starting point:
//
//      X = [ -0.5, 0.25 ] (easy convergence)
//      X = [ -4.0, 5.00 ] (harder convergence)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Zbigniew Michalewicz,
//    Genetic Algorithms + Data Structures = Evolution Programs,
//    Third Edition,
//    Springer Verlag, 1996,
//    ISBN: 3-540-60676-9,
//    LC: QA76.618.M53.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double GOLDSTEIN_PRICE, the value of the function at X.
//
{
  double a;
  double b;
  double c;
  double d;
  double f;

  a = x[0] + x[1] + 1.0;

  b = 19.0 - 14.0 * x[0] + 3.0 * x[0] * x[0] - 14.0 * x[1] 
    + 6.0 * x[0] * x[1] + 3.0 * x[1] * x[1];

  c = 2.0 * x[0] - 3.0 * x[1];

  d = 18.0 - 32.0 * x[0] + 12.0 * x[0] * x[0] + 48.0 * x[1] 
    - 36.0 * x[0] * x[1] + 27.0 * x[1] * x[1];

  f = ( 1.0 + a * a * b ) * ( 30.0 + c * c * d );

  return f;
}
//****************************************************************************80

double himmelblau ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HIMMELBLAU computes the Himmelblau function.
//
//  Discussion:
//
//    This function has 4 global minima:
//
//      X* = (  3,        2       ), F(X*) = 0.
//      X* = (  3.58439, -1.84813 ), F(X*) = 0.
//      X* = ( -3.77934, -3.28317 ), F(X*) = 0.
//      X* = ( -2.80512,  3.13134 ), F(X*) = 0.
//
//    Suggested starting points are
//
//      (+1,+1),
//      (-1,+1),
//      (-1,-1),
//      (+1,-1),
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Himmelblau,
//    Applied Nonlinear Programming,
//    McGraw Hill, 1972,
//    ISBN13: 978-0070289215,
//    LC: T57.8.H55.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double HIMMELBLAU, the value of the function at X.
//
{
  double f;

  f = pow ( x[0] * x[0] + x[1] - 11.0, 2 ) 
    + pow ( x[0] + x[1] * x[1] - 7.0, 2 );
 
  return f;
}
//****************************************************************************80

double local ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LOCAL computes the local function.
//
//  Discussion:
//
//    This function has a local minimizer:
//
//      X* = ( 0.2858054412&, 0.2793263206&), F(X*) = 5.9225&
//
//    and a global minimizer:
//
//      X* = ( -21.02667179&, -36.75997872&), F(X*) = 0.
//
//    Suggested starting point for local minimizer:
//
//      X = ( 1, 1 ), F(X) = 3.33 * 10^6.
//
//    Suggested starting point for global minimizer:
//
//      X = ( -15, -35), F(X) = 1.49 * 10^8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Himmelblau,
//    Applied Nonlinear Programming,
//    McGraw Hill, 1972,
//    ISBN13: 978-0070289215,
//    LC: T57.8.H55.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double LOCAL, the value of the function at X.
//
{
  double f;

  f = pow ( x[0] * x[0] + 12.0 * x[1] - 1.0, 2 )
    + pow ( 49.0 * x[0] * x[0] + 49.0 * x[1] * x[1] + 84.0 * x[0] 
    + 2324.0 * x[1] - 681.0, 2 );
 
  return f;
}
//****************************************************************************80

double mckinnon ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MCKINNON computes the McKinnon function.
//
//  Discussion:
//
//    This function has a global minimizer:
//
//      X* = ( 0.0, -0.5 ), F(X*) = -0.25.
//
//    There are three parameters, TAU, THETA and PHI.
//
//    1 < TAU, then F is strictly convex.
//             and F has continuous first derivatives.
//    2 < TAU, then F has continuous second derivatives.
//    3 < TAU, then F has continuous third derivatives.
//
//    However, this function can cause the Nelder-Mead optimization
//    algorithm to "converge" to a point which is not the minimizer
//    of the function F.
//
//    Sample parameter values which cause problems for Nelder-Mead 
//    include:
//
//      PHI = 10,  TAU = 1, THETA = 15
//      PHI = 60,  TAU = 2, THETA =  6
//      PHI = 400, TAU = 3, THETA =  6
//
//    To get the bad behavior, we also assume the initial simplex has the form
//
//      X1 = (0,0),
//      X2 = (1,1),
//      X3 = (A,B), 
//
//    where 
//
//      A = (1+sqrt(33))/8 =  0.84307&
//      B = (1-sqrt(33))/8 = -0.59307&
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ken McKinnon,
//    Convergence of the Nelder-Mead simplex method to a nonstationary point,
//    SIAM Journal on Optimization,
//    Volume 9, Number 1, 1998, pages 148-158.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double MCKINNON, the value of the function at X.
//
{
  double f;
  double phi;
  double tau;
  double theta;

  mckinnon_parameters ( "get", phi, tau, theta );

  if ( x[0] <= 0.0 )
  {
    f = theta * phi * pow ( r8_abs ( x[0] ), tau ) + x[1] * ( 1.0 + x[1] );
  }
  else
  {
    f = theta       * pow (           x[0],  tau ) + x[1] * ( 1.0 + x[1] );
  }

  return f;
}
//****************************************************************************80

void mckinnon_parameters ( string action, double &phi, double &tau, double &theta )

//****************************************************************************80
//
//  Purpose:
//
//    MCKINNON_PARAMETERS sets or gets McKinnon function parameters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ACTION.  
//    "S", set the parameters.
//    "G", get the parameters.
//
//    Input/output, double &PHI, &TAU, &THETA, the parameter values.
//
{
  static double phi_save = 60.0;
  static double tau_save = 2.0;
  static double theta_save = 6.0;

  if ( action[0] == 'S' || action[0] == 's' ) 
  {
    phi_save = phi;
    tau_save = tau;
    theta_save = theta;
  }
  else if ( action[0] == 'G' || action[0] == 'g' )
  {
    phi = phi_save;
    tau = tau_save;
    theta = theta_save;
  }
  else
  {
    cout << "\n";
    cout << "MCKINNON_PARAMETERS - Fatal error!\n";
    cout << "  Unexpected value of ACTION.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

double powell ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    POWELL computes the Powell singular quartic function.
//
//  Discussion:
//
//    This function has a global minimizer:
//
//      X* = ( 0.0, 0.0, 0.0, 0.0 ), F(X*) = 0.
//
//    Start the search at
//
//      X = ( 3.0, -1.0, 0.0, 1.0 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Michael Powell,
//    An Iterative Method for Finding Stationary Values of a Function
//    of Several Variables,
//    Computer Journal,
//    Volume 5, 1962, pages 147-151.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double POWELL, the value of the function at X.
//
{
  double f;
  double f1;
  double f2;
  double f3;
  double f4;

  f1 = x[0] + 10.0 * x[1];
  f2 =                            x[2] - x[3];
  f3 =               x[1] - 2.0 * x[2];
  f4 = x[0]                            - x[3];

  f = f1 * f1 + f2 * f2 + f3 * f3 + f4 * f4;
 
  return f;
}
//****************************************************************************80

double rosenbrock ( int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK computes the Rosenbrock function.
//
//  Discussion:
//
//    There is a global minimum at X* = (1,1), F(X*) = 0.
//
//    The starting point X = [ -1.2, 1.0 ] is suggested.
//
//    The contours are sharply twisted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Rosenbrock,
//    An Automatic Method for Finding the Greatest or Least Value of a Function,
//    Computer Journal,
//    Volume 3, 1960, pages 175-184.
//
//  Parameters:
//
//    Input, int M, the number of variables.
//
//    Input, double X[M], the argument of the function.
//
//    Output, double ROSENBROCK, the value of the function at X.
//
{
  double f;

  f = pow ( 1.0 - x[0], 2 ) + 100.0 * pow ( x[1] - x[0] * x[0], 2 );

  return f;
}
