# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "lagrange_approx_1d.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"
# include "test_interp_1d.hpp"

int main ( );
void test02 ( int prob, int m, int nd );
void test03 ( int prob, int m, int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGRANGE_APPROX_1D_PRB.
//
//  Discussion:
//
//    LAGRANGE_APPROX_1D_PRB tests the LAGRANGE_APPROX_1D library.
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
  int m;
  int m_test[7] = { 0, 1, 2, 3, 4, 8, 16 };
  int m_test_num = 7;
  int nd;
  int nd_test[3] = { 16, 64, 1000 };
  int nd_test_num = 3;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "LAGRANGE_APPROX_1D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGRANGE_APPROX_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  These tests need the TEST_INTERP_1D library.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      for ( k = 0; k < nd_test_num; k++ )
      {
        nd = nd_test[k];
        test02 ( prob, m, nd );
      }
    }
  }

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      for ( k = 0; k < nd_test_num; k++ )
      {
        nd = nd_test[k];
        test03 ( prob, m, nd );
      }
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGRANGE_APPROX_1D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test02 ( int prob, int m, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests LAGRANGE_APPROX_1D with evenly spaced data
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
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int M, the polynomial approximant degree.
//
//    Input, int ND, the number of data points.
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
  cout << "TEST02:\n";
  cout << "  Approximate evenly spaced data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use polynomial approximant of degree " << m << "\n";
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
//  #1:  Does approximant come close to function at data points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = lagrange_approx_1d ( m, nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per data node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
//****************************************************************************80

void test03 ( int prob, int m, int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests LAGRANGE_APPROX_1D with Chebyshev spaced data.
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
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int M, the polynomial approximant degree.
//
//    Input, int ND, the number of data points.
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
  cout << "TEST03:\n";
  cout << "  Approximate Chebyshev-spaced data from TEST_INTERP_1D problem #" << prob << "\n";
  cout << "  Use polynomial approximant of degree " << m << "\n";
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
  yi = lagrange_approx_1d ( m, nd, xd, yd, ni, xi );

  int_error = r8vec_norm_affine ( nd, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per data node = " << int_error << "\n";

  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;

  return;
}
