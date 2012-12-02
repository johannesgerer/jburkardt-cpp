# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "shepard_interp_1d.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, double p );

//****************************************************************************80

int main ( )

/******************************************************************************/
//
//  Purpose:
//
//    SHEPARD_INTERP_1D_TEST tests SHEPARD_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  double p;
  int p_num = 5;
  double p_test[5] = { 0.0, 1.0, 2.0, 4.0, 8.0 };
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "SHEPARD_INTERP_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the SHEPARD_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test needs the TEST_INTERP library as well.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < p_num; j++ )
    {
      p = p_test[j];
      test01 ( prob, p );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "SHEPARD_INTERP_1D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, double p )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SHEPARD_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  double int_error;
  int i;
  double ld;
  double li;
  int nd;
  int ni;
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
  cout << "  using Shepard interpolation with P = " << p << "\n";

  dim_num = p00_dim_num ( prob );

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, dim_num, nd );
  
  if ( p == 0.0 )
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
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = new double[ni];
  for ( i = 0; i < ni; i++ )
  {
    xi[i] = xd[i];
  }

  yi = shepard_interp_1d ( nd, xd, yd, p, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " 
    << int_error << "\n";

  delete [] xi;
  delete [] yi;
//
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
//
  xmin = r8vec_min ( nd, xd );
  xmax = r8vec_max ( nd, xd );
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, xmin, xmax );
  yi = shepard_interp_1d ( nd, xd, yd, p, ni, xi );

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

  cout << "\n";
  cout << "  Normalized length of piecewise linear interpolant = " << ld << "\n";
  cout << "  Normalized length of Shepard interpolant          = " << li << "\n";

  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
