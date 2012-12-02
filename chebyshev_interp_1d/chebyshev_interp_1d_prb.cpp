# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

using namespace std;

# include "chebyshev_interp_1d.hpp"
# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_INTERP_1D_TEST tests CHEBYSHEV_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2012
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
  cout << "CHEBYSHEV_INTERP_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the CHEBYSHEV_INTERP_1D library.\n";
  cout << "  The QR_SOLVE and R8LIB libraries are needed.\n";
  cout << "  The test needs the TEST_INTERP library.\n";

  prob_num = p00_prob_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test01 ( prob );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "CHEBYSHEV_INTERP_1D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CHEBYSHEV_VALUE_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2012
//
//  Author:
//
//    John Burkardt
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
  cout << "CHEBYSHEV_INTERP_1D_TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP problem #" << prob << "\n";

  nd = p00_data_num ( prob );
  cout << "  Number of data points = " << nd << "\n";

  xy = p00_data ( prob, 2, nd );
  
  r8mat_transpose_print ( 2, nd, xy, "  Data array:" );

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
  yi = chebyshev_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
