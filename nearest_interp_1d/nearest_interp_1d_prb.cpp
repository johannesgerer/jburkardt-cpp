# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "nearest_interp_1d.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void nearest_interp_1d_test01 ( int prob, int ni );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    NEAREST_INTERP_1D_TEST tests NEAREST_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2012
//
//  Author:
//
//   John Burkardt
//
{
  int ni;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "NEAREST_INTERP_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the NEAREST_INTERP_1D library.\n";
  cout << "  The test needs the TEST_INTERP library.\n";

  prob_num = p00_prob_num ( );

  ni = 11;
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    nearest_interp_1d_test01 ( prob, ni );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "NEAREST_INTERP_1D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void nearest_interp_1d_test01 ( int prob, int ni )

//****************************************************************************80
//
//  Purpose:
//
//    NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the problem.
//
//    Input, int NI, the number of interpolation points.
//
{
  double *d;
  int j;
  int nd;
  char title_ch[80];
  string title;
  double *xd;
  double *xi;
  double xd_max;
  double xd_min;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "NEAREST_INTERP_1D_TEST01\n";
  cout << "  Sample the nearest neighbor interpolant for problem # " << prob << "\n";

  nd = p00_data_num ( prob );

  d = p00_data ( prob, 2, nd );

  xd = new double[nd];
  yd = new double[nd];

  for ( j = 0; j < nd; j++ )
  {
    xd[j] = d[0+j*2];
    yd[j] = d[1+j*2];
  }

  xd_min = r8vec_min ( nd, xd );
  xd_max = r8vec_max ( nd, xd );

  xi = r8vec_linspace_new ( ni, xd_min, xd_max );
  yi = nearest_interp_1d ( nd, xd, yd, ni, xi );

  sprintf ( title_ch, "X, Y for problem %d", prob );
  title = title_ch;

  r8vec2_print ( ni, xi, yi, title );

  delete [] d;
  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}

