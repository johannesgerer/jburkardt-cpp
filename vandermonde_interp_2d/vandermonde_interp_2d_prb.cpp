# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "vandermonde_interp_2d.hpp"
# include "qr_solve.hpp"
# include "test_interp_2d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int m );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_INTERP_2D_PRB tests VANDERMONDE_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int m;
  int m_test[5] = { 1, 2, 3, 4, 8 };
  int m_test_num = 5;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "VANDERMONDE_INTERP_2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the VANDERMONDE_INTERP_2D library.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      test01 ( prob, m );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "VANDERMONDE_INTERP_2D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int m )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests VANDERMONDE_INTERP_2D_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the degree of interpolation.
//
{
  double *a;
  double app_error;
  double *c;
  bool debug = false;
  int nd;
  int ni;
  int seed;
  int tmp1;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double *zd;
  double *zi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_2D problem #" << prob << "\n";
  cout << "  Create an interpolant of total degree " << m << "\n";
  tmp1 = triangle_num ( m + 1 );
  cout << "  Number of data values needed is " << tmp1 << "\n";

  nd = tmp1;

  seed = 123456789;

  xd = r8vec_uniform_01_new ( nd, seed );
  yd = r8vec_uniform_01_new ( nd, seed );
  zd = new double[nd];
  f00_f0 ( prob, nd, xd, yd, zd );

  if ( debug )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
//
//  Compute the Vandermonde matrix.
//
  a = vandermonde_interp_2d_matrix ( nd, m, xd, yd );
//
//  Solve linear system.
//
  c = qr_solve ( nd, nd, a, zd );
//
//  #1:  Does interpolant match function at data points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );
  zi = r8poly_value_2d ( m, c, ni, xi, yi );

  app_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 data interpolation error = " << app_error << "\n";

  delete [] a;
  delete [] c;
  delete [] xd;
  delete [] xi;
  delete [] yd;
  delete [] yi;
  delete [] zd;
  delete [] zi;

  return;
}
