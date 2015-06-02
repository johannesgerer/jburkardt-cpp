# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "vandermonde_approx_2d.hpp"
# include "test_interp_2d.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int grid, int m );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for VANDERMONDE_APPROX_2D_PRB.
//
//  Discussion:
//
//    VANDERMONDE_APPROX_2D_PRB tests the VANDERMONDE_APPROX_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int m;
  int m_test[5] = { 0, 1, 2, 4, 8 };
  int m_test_num = 5;
  int grid;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "VANDERMONDE_APPROX_2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the VANDERMONDE_APPROX_2D library.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test also needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    grid = 1;
    for ( j = 0; j < m_test_num; j++ )
    {
      m = m_test[j];
      test01 ( prob, grid, m );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "VANDERMONDE_APPROX_2D_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int grd, int m )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_2D_TEST01 tests VANDERMONDE_APPROX_2D_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int GRD, the grid number.
//    (Can't use GRID as the name because that's also a plotting function.)
//
//    Input, int M, the total polynomial degree.
//
{
  double *a;
  double app_error;
  double *c;
  int nd;
  int ni;
  int tm;
  double *xd;
  double *xi;
  double *yd;
  double *yi;
  double *zd;
  double *zi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Approximate data from TEST_INTERP_2D problem #" << prob << "\n";
  cout << "  Use grid from TEST_INTERP_2D with index #" << grd << "\n";
  cout << "  Using polynomial approximant of total degree " << m << "\n";

  nd = g00_size ( grd );
  cout << "  Number of data points = " << nd << "\n";

  xd = new double[nd];
  yd = new double[nd];
  g00_xy ( grd, nd, xd, yd );

  zd = new double[nd];
  f00_f0 ( prob, nd, xd, yd, zd );

  if ( nd < 10 )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
//
//  Compute the Vandermonde matrix.
//
  tm = triangle_num ( m + 1 );
  a = vandermonde_approx_2d_matrix ( nd, m, tm, xd, yd );
//
//  Solve linear system.
//
  c = qr_solve ( nd, tm, a, zd );
//
//  #1:  Does approximant match function at data points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );
  zi = r8poly_value_2d ( m, c, ni, xi, yi );

  app_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 data approximation error = " << app_error << "\n";

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
