# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "shepard_interp_nd.hpp"
# include "test_interp_nd.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, double p, int m, int nd );
void test02 ( int prob, double p, int m, int n1d );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SHEPARD_INTERP_ND_PRB.
//
//  Discussion:
//
//    SHEPARD_INTERP_ND_PRB tests the SHEPARD_INTERP_ND library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int m;
  int n1d;
  int nd;
  double p;
  double p_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  int p_test_num = 4;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "SHEPARD_INTERP_ND_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the SHEPARD_INTERP_ND library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test also needs the TEST_INTERP_ND library.\n";
//
//  Look at Shepard interpolant on an irregular grid.
//
  nd = 25;

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( m = 2; m <= 5; m = m + 3 )
    {
      for ( j = 0; j < p_test_num; j++ )
      {
        p = p_test[j];
        test01 ( prob, p, m, nd );
      }

    }
  }
//
//  Look at Shepard interpolant on a regular N1D^M grid.
//
  n1d = 5;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( m = 2; m <= 5; m = m + 3 )
    {
      for ( j = 0; j < p_test_num; j++ )
      {
        p = p_test[j];
        test02 ( prob, p, m, n1d );
      }
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "SHEPARD_INTERP_ND_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, double p, int m, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SHEPARD_INTERP on an irregular grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, double P, the power used in the distance weighting.
//
//    Input, int M, the spatial dimension.
//
//    Input, int ND, the number of data points.
//
{
  double app_error;
  double *c;
  int i;
  double int_error;
  int j;
  int ni;
  int seed;
  double *w;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_ND problem #" << prob << "\n";
  cout << "  using Shepard interpolation with P = " << p  << "\n";
  cout << "  spatial dimension M = " << m  << "\n";
  cout << "  and an irregular grid of ND = " << nd << " data points.\n";
//
//  Set problem parameters:
//
  seed = 123456789;
  c = r8vec_uniform_01_new ( m, seed );
  w = r8vec_uniform_01_new ( m, seed );

  xd = r8mat_uniform_01_new ( m, nd, seed );

  zd = p00_f ( prob, m, c, w, nd, xd );
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] zi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  ni = 50;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );
  ze = p00_f ( prob, m, c, w, ni, xi );

  app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

  cout << "  L2 approximation error averaged per 1000 samples =     " << app_error << "\n";

  delete [] c;
  delete [] w;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

void test02 ( int prob, double p, int m, int n1d )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests SHEPARD_INTERP_ND on a regular N1D^M grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, double P, the power used in the distance weighting.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N1D, the number of points in 1D.
//
{
  double a;
  double app_error;
  double b;
  double *c;
  int i;
  double int_error;
  int nd;
  int ni;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;
//
//  Set problem parameters:
//
  seed = 123456789;
  c = r8vec_uniform_01_new ( m, seed );
  w = r8vec_uniform_01_new ( m, seed );

  nd = i4_power ( n1d, m );

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Interpolate data from TEST_INTERP_ND problem #" << prob << "\n";
  cout << "  using Shepard interpolation with P = " << p << "\n";
  cout << "  spatial dimension M = " << m << "\n";
  cout << "  and a regular grid of N1D^M = " << nd << " data points.\n";

  a = 0.0;
  b = 1.0;

  x1d = r8vec_linspace_new ( n1d, a, b );

  xd = new double[m*nd];
  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  zd = p00_f ( prob, m, c, w, nd, xd );
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, nd, xd );
  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] zi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  xi = r8mat_uniform_01_new ( m, ni, seed );

  zi = shepard_interp_nd ( m, nd, xd, zd, p, ni, xi );

  ze = p00_f ( prob, m, c, w, ni, xi );

  app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

  cout << "  L2 approximation error averaged per 1000 samples =     " << app_error << "\n";

  delete [] c;
  delete [] w;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
