# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "vandermonde_approx_1d.hpp"
# include "test_interp.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int m );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_1D_TEST tests VANDERMONDE_APPROX_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int m;
  int m_test[8] = { 0, 1, 2, 3, 4, 5, 9, 12 };
  int m_test_num = 8;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "VANDERMONDE_APPROX_1D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the VANDERMONDE_APPROX_1D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The QR_SOLVE library is needed.\n";
  cout << "  The test needs the CONDITION library.\n";
  cout << "  The test needs the TEST_INTERP libary.\n";

  prob_num = p00_prob_num ( );
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
  cout << "VANDERMONDE_APPROX_1D_TEST:\n";
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
//    TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the polynomial degree.
//
{
  double *a;
  double app_error;
  double *c;
  int debug = 0;
  int i;
  int j;
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
  cout << "  Approximate data from TEST_INTERP problem #" << prob << "\n";

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
//
//  Compute the Vandermonde matrix.
//
  cout << "  Using polynomial approximant of degree " << m << "\n";

  a = vandermonde_approx_1d_matrix ( nd, m, xd );
//
//  Solve linear system.
//
  c = qr_solve ( nd, m + 1, a, yd );
//
//  #1:  Does approximant match function at data points?
//
  ni = nd;
  xi = r8vec_copy_new ( ni, xd );
  yi = r8poly_value ( m, c, ni, xi );

  app_error = r8vec_norm_affine ( ni, yi, yd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 data approximation error = " << app_error << "\n";

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

  yi = r8poly_value ( m, c, ni, xi );

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
  cout << "  Normalized length of polynomial interpolant       = " << li << "\n";

  delete [] a;
  delete [] c;
  delete [] xd;
  delete [] xi;
  delete [] xy;
  delete [] yd;
  delete [] yi;

  return;
}
