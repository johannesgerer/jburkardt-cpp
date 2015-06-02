# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "rbf_interp_1d.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, void phi ( int n, double r[], double r0, double v[] ), 
  string phi_name, double r0 );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RBF_INTERP_1D_PRB.
//
//  Discussion:
//
//    RBF_INTERP_1D_PRB tests the RBF_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nd;
  int prob;
  int prob_num;
  double r0;
  double *xd;
  double xmax;
  double xmin;
  double *xy;

  timestamp ( );
  cout << "\n";
  cout << "RBF_INTERP_1D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the RBF_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The test needs the TEST_INTERP library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Determine an appropriate value of R0, the spacing parameter.
//
    nd = p00_data_num ( prob );
    xy = p00_data ( prob, 2, nd );
    xd = ( double * ) malloc ( nd * sizeof ( double ) );
    for ( i = 0; i < nd; i++ )
    {
      xd[i] = xy[0+i*2];
    }
    xmax = r8vec_max ( nd, xd );
    xmin = r8vec_min ( nd, xd );
    r0 = ( xmax - xmin ) / ( double ) ( nd - 1 );
    delete [] xd;
    delete [] xy;
    test01 ( prob, phi1, "phi1", r0 );
    test01 ( prob, phi2, "phi2", r0 );
    test01 ( prob, phi3, "phi3", r0 );
    test01 ( prob, phi4, "phi4", r0 );
  }
/*
  Terminate.
*/
  cout << "\n";
  cout << "RBF_INTERP_1D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, void phi ( int n, double r[], double r0, double v[] ), 
  string phi_name, double r0 )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests RBF_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the problem.
//
//    Input, double PHI ( int n, double r[], double r0, double v[] ), 
//    the name of the radial basis function.
//
//    Input, string PHI_NAME, the name of the radial basis function.
//
//    Input, double R0, the scale factor.  Typically, this might be
//    a small multiple of the average distance between points.
//
{
  bool debug = false;
  int i;
  double int_error;
  double ld;
  double li;
  int m;
  int nd;
  int ni;
  double *w;
  double *xd;
  double *xi;
  double xmax;
  double xmin;
  double *xy;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP problem #" << prob << "\n";
  cout << "  using radial basis function \"" << phi_name << "\".\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, 2, nd );
  
  if ( debug )
  {
    r8mat_transpose_print ( 2, nd, xy, "  Data array:" );
  }

  xd = new double[nd];
  yd = new double[nd];
  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+i*2];
    yd[i] = xy[1+i*2];
  }
  m = 1;
  w = rbf_weight ( m, nd, xd, r0, phi, yd );
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = rbf_interp ( m, nd, xd, r0, phi, w, ni, xi );
  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] yi;
//
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
//
  xmax = r8vec_max ( nd, xd );
  xmin = r8vec_min ( nd, xd );
  ymax = r8vec_max ( nd, yd );
  ymin = r8vec_min ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = rbf_interp ( m, nd, xd, r0, phi, w, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) );
  }
  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( xmax - xmin ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  cout << "  Normalized length of piecewise linear interpolant = " << ld << "\n";
  cout << "  Normalized length of polynomial interpolant       = " << li << "\n";

  delete [] w;
  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
