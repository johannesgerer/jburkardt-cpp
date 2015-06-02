# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "barycentric_interp_1d.hpp"
# include "test_interp_1d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int nd );
void test02 ( int prob, int nd );
void test03 ( int prob, int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BARYCENTRIC_INTERP_1D_PRB.
//
//  Discussion:
//
//    BARYCENTRIC_INTERP_1D_PRB tests the BARYCENTRIC_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int nd;
  int nd_test[6] = { 4, 8, 16, 32, 64, 1000 };
  int nd_test_num = 6;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "BARYCENTRIC_INTERP_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the BARYCENTRIC_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The tests need the TEST_INTERP_1D library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test01 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test02 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < nd_test_num; i++ )
    {
      nd = nd_test[i];
      test03 ( prob, nd );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "BARYCENTRIC_INTERP_1D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    BARYCENTRIC_INTERP_1D_TEST01 tests LAGCHEBY1_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int ND, the number of data points to use.
//
{
  double a;
  double b;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "BARYCENTRIC_INTERP_1D_TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use Chebyshev Type 1 spacing for data points.\n";
  cout << "  Number of data points = " << nd << "\n";
//
//  Define the data.
//
  a =  0.0;
  b = +1.0;
  xd = r8vec_cheby1space_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  #1:  Does the interpolant match the function at the interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagcheby1_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

void test02 ( int prob, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    BARYCENTRIC_INTERP_1D_TEST02 tests LAGCHEBY2_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int ND, the number of data points to use.
//
{
  double a;
  double b;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "BARYCENTRIC_INTERP_1D_TEST02:\n";
  cout << "  Interpolate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use Chebyshev Type 2 spacing for data points.\n";
  cout << "  Number of data points = " << nd << "\n";
//
//  Define the data.
//
  a =  0.0;
  b = +1.0;
  xd = r8vec_cheby2space_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  #1:  Does the interpolant match the function at the interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagcheby2_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

void test03 ( int prob, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    BARYCENTRIC_INTERP_1D_TEST03 tests LAGEVEN_INTERP_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int ND, the number of data points to use.
//
{
  double a;
  double b;
  int i;
  double int_error;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;

  cout << "\n";
  cout << "BARYCENTRIC_INTERP_1D_TEST03:\n";
  cout << "  Interpolate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use even spacing for data points.\n";
  cout << "  Number of data points = " << nd << "\n";
//
//  Define the data.
//
  a =  0.0;
  b = +1.0;
  xd = r8vec_midspace_new ( nd, a, b );
  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  #1:  Does the interpolant match the function at the interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lageven_interp_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
