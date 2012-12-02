# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_1d.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void pwl_interp_1d_test01 ( int prob );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_1D_TEST tests PWL_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "PWL_INTERP_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the PWL_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The test needs the TEST_INTERP library.\n";

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    pwl_interp_1d_test01 ( prob );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PWL_INTERP_1D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void pwl_interp_1d_test01 ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_1D_TEST01 tests PWL_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
{
  int i;
  double int_error;
  int nd;
  int ni;
  double *xd;
  double *xi;
  double *xy;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "PWL_INTERP_1D_TEST01:\n";
  cout << "  PWL_INTERP_1D evaluates the piecewise linear interpolant.\n";
  cout << "  Interpolate data from TEST_INTERP problem #" << prob << "\n";

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, 2, nd );
  
  r8mat_transpose_print ( 2, nd, xy, "  Data array:" );

  xd = new double[nd];
  yd = new double[nd];

  for ( i = 0; i < nd; i++ )
  {
    xd[i] = xy[0+2*i];
    yd[i] = xy[1+2*i];
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

  yi = pwl_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_diff_norm ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
