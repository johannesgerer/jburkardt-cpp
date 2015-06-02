# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_1d.hpp"
# include "test_interp_1d.hpp"
# include "r8lib.hpp"

int main ( );
void test02 ( int prob, int nd );
void test03 ( int prob, int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGRANGE_INTERP_1D_PRB.
//
//  Discussion:
//
//    LAGRANGE_INTERP_1D_PRB tests the LAGRANGE_INTERP_1D library.
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
  int nd_test_num = 6;

  int j;
  int nd;
  int nd_test[6] = { 4, 8, 16, 32, 64, 256 };
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "LAGRANGE_INTERP_1D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGRANGE_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  These tests need the TEST_INTERP_1D library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < nd_test_num; j++ )
    {
      nd = nd_test[j];
      test02 ( prob, nd );
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < nd_test_num; j++ )
    {
      nd = nd_test[j];
      test03 ( prob, nd );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGRANGE_INTERP_1D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( int prob, int nd )

//*****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests LAGRANGE_VALUE_1D with evenly spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2012
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
  double ld;
  double li;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Interpolate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use even spacing for data points.\n";
  cout << "  Number of data points = " << nd << "\n";

  a = 0.0;
  b = 1.0;
  
  xd = r8vec_linspace_new ( nd, a, b );

  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] yi;
//
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
//
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, a, b );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( b - a ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( b - a ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  cout << "\n";
  cout << "  Normalized length of piecewise linear interpolant = " << ld << "\n";
  cout << "  Normalized length of polynomial interpolant       = " << li << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

void test03 ( int prob, int nd )

//*****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests LAGRANGE_VALUE_1D with Chebyshev spaced data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2012
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
  double ld;
  double li;
  int ni;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double ymax;
  double ymin;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Interpolate data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use Chebyshev spacing for data points.\n";
  cout << "  Number of data points = " << nd << "\n";

  a = 0.0;
  b = 1.0;
  
  xd = r8vec_chebyspace_new ( nd, a, b );

  yd = p00_f ( prob, nd, xd );

  if ( nd < 10 )
  {
    r8vec2_print ( nd, xd, yd, "  Data array:" );
  }
//
//  #1:  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] xi;
  delete [] yi;
//
//  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
//  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
//  (YMAX-YMIN).
//
  ymin = r8vec_min ( nd, yd );
  ymax = r8vec_max ( nd, yd );

  ni = 501;
  xi = r8vec_linspace_new ( ni, a, b );
  yi = lagrange_value_1d ( nd, xd, yd, ni, xi );

  ld = 0.0;
  for ( i = 0; i < nd - 1; i++ )
  {
    ld = ld + sqrt ( pow ( ( xd[i+1] - xd[i] ) / ( b - a ), 2 )
                   + pow ( ( yd[i+1] - yd[i] ) / ( ymax - ymin ), 2 ) ); 
  }

  li = 0.0;
  for ( i = 0; i < ni - 1; i++ )
  {
    li = li + sqrt ( pow ( ( xi[i+1] - xi[i] ) / ( b - a ), 2 )
                   + pow ( ( yi[i+1] - yi[i] ) / ( ymax - ymin ), 2 ) );
  }

  cout << "\n";
  cout << "  Normalized length of piecewise linear interpolant = " << ld << "\n";
  cout << "  Normalized length of polynomial interpolant       = " << li << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
