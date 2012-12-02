# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "shepard_interp_2d.hpp"
# include "test_interp_2d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int g, double p );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_INTERP_2D_TEST tests SHEPARD_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int g;
  int j;
  double p;
  double p_test[4] = { 1.0, 2.0, 4.0, 8.0 };
  int p_test_num = 4;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "SHEPARD_INTERP_2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the SHEPARD_INTERP_2D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test also needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );
  g = 1;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < p_test_num; j++ )
    {
      p = p_test[j];
      test01 ( prob, g, p );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "SHEPARD_INTERP_2D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int g, double p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SHEPARD_INTERP_2D.
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
//    Input, int G, the grid number.
//
//    Input, double P, the power used in the distance weighting.
//
{
  bool debug = false;
  double int_error;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double *zd;
  double *zi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_2D problem #" << prob << "\n";
  cout << "  using grid #" << g << "\n";
  cout << "  using Shepard interpolation with P = " << p << "\n";

  nd = g00_size ( g );
  cout << "  Number of data points = " << nd << "\n";

  xd = new double[nd];
  yd = new double[nd];
  g00_xy ( g, nd, xd, yd );
  
  zd = new double[nd];
  f00_f0 ( prob, nd, xd, yd, zd );

  if ( debug )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );

  zi = shepard_interp_2d ( nd, xd, yd, zd, p, ni, xi, yi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;
  delete [] zd;
  delete [] zi;

  return;
}
