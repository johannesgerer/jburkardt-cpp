# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "rbf_interp_2d.hpp"
# include "r8lib.hpp"
# include "test_interp_2d.hpp"

int main ( );
void test01 ( int prob, int g, 
  void phi ( int n, double r[], double r0, double v[] ), string phi_name );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RBF_INTERP_2D_PRB.
//
//  Discussion:
//
//    RBF_INTERP_2D_PRB tests the RBF_INTERP_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int g;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "RBF_INTERP_2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the RBF_INTERP_2D library.\n";
  cout << "  The R8LIB library is required.\n";
  cout << "  This test also needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );
  g = 1;

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob, g, phi1, "phi1" );
    test01 ( prob, g, phi2, "phi2" );
    test01 ( prob, g, phi3, "phi3" );
    test01 ( prob, g, phi4, "phi4" );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "RBF_INTERP_2D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int g, 
  void phi ( int n, double r[], double r0, double v[] ), string phi_name )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_2D_TEST01 tests RBF_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the problem.
//
//    Input, int G, the index of the grid.
//
//    Input, void PHI ( int n, double r[], double r0, double v[] ), the 
//    radial basis function.
//
//    Input, string PHI_NAME, the name of the radial basis function.
//
{
  int debug = 0;
  double e;
  int i;
  double int_error;
  int m;
  int nd;
  int ni;
  double r0;
  double volume;
  double *w;
  double *xd;
  double xmax;
  double xmin;
  double *xyd;
  double *xyi;
  double *yd;
  double ymax;
  double ymin;
  double *zd;
  double *zi;

  cout << "\n";
  cout << "RBF_INTERP_2D_TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_2D problem #" << prob << "\n";
  cout << "  using grid #" << g << "\n";
  cout << "  using radial basis function \"" << phi_name << "\".\n";

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

  m = 2;
  xyd = new double[2*nd];

  for ( i = 0; i < nd; i++ )
  {
    xyd[0+i*2] = xd[i];
    xyd[1+i*2] = yd[i];
  }

  xmax = r8vec_max ( nd, xd );
  xmin = r8vec_min ( nd, xd );
  ymax = r8vec_max ( nd, yd );
  ymin = r8vec_min ( nd, yd );
  volume = ( xmax - xmin ) * ( ymax - ymin );

  e = 1.0 / ( double ) ( m );
  r0 = pow ( volume / nd, e );

  cout << "  Setting R0 = " << r0 << "\n";

  w = rbf_weight ( m, nd, xyd, r0, phi, zd );
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xyi = r8mat_copy_new ( 2, ni, xyd );

  zi = rbf_interp ( m, nd, xyd, r0, phi, w, ni, xyi );

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " 
       << int_error << "\n";

  delete [] w;
  delete [] xd;
  delete [] xyd;
  delete [] xyi;
  delete [] yd;
  delete [] zd;
  delete [] zi;

  return;
}
