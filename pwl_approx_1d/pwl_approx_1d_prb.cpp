# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_interp_1d.hpp"
# include "pwl_approx_1d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int nc, int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PWL_APPROX_1D_PRB.
//
//  Discussion:
//
//    PWL_APPROX_1D_PRB tests the PWL_APPROX_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int k;
  int nc;
  int nc_test[4] = { 2, 4, 8, 16 };
  int nc_test_num = 4;
  int nd;
  int nd_test[2] = { 16, 64 };
  int nd_test_num = 2;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "PWL_APPROX_1D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PWL_APPROX_1D library.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The test also needs the TEST_INTERP_1D library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < nc_test_num; j++ )
    {
      nc = nc_test[j];
      for ( k = 0; k < nd_test_num; k++ )
      {
        nd = nd_test[k];
        test01 ( prob, nc, nd );
      }
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PWL_APPROX_1D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int nc, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PWL_APPROX_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int NC, the number of control points.
//
//    Input, int ND, the number of data points.
//
{
  double app_error;
  int ni;
  double *xc;
  double *xd;
  double *xi;
  double *yc;
  double *yd;
  double *yi;
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Approximate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Number of control points = " << nc << "\n";
  cout << "  Number of data points = " << nd << "\n";

  xmin = 0.0;
  xmax = 1.0;

  xd = r8vec_linspace_new ( nd, xmin, xmax );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  Determine control values.
// 
  xc = r8vec_linspace_new ( nc, xmin, xmax );
  yc = pwl_approx_1d ( nd, xd, yd, nc, xc );
//
//  #1:  Does approximant come close to function at data points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = pwl_interp_1d ( nc, xc, yc, ni, xi );

  app_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per data node = " << app_error << "\n";

  delete [] xc;
  delete [] xd;
  delete [] xi;
  delete [] yc;
  delete [] yd;
  delete [] yi;

  return;
}
