# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_2d.hpp"
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
//    MAIN is the main program for LAGRANGE_INTERP_2D_PRB.
//
//  Discussion:
//
//    LAGRANGE_INTERP_2D_PRB tests the LAGRANGE_INTERP_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  int m_test[5] = { 1, 2, 3, 4, 8 };
  int m_test_num = 5;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "LAGRANGE_INTERP_2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGRANGE_INTERP_2D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test also needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );
//
//  Numerical tests.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < m_test_num; i++ )
    {
      m = m_test[i];
      test01 ( prob, m );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGRANGE_INTERP_2D_PRB:\n";
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
//    LAGRANGE_INTERP_2D_TEST01 tests LAGRANGE_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the polynomial degree in each dimension.
//
{
  double app_error;
  int i;
  int ij;
  double int_error;
  int j;
  int mx;
  int my;
  int nd;
  int ni;
  double *xd;
  double *xd_1d;
  double *xi;
  double *xi_1d;
  double *yd;
  double *yd_1d;
  double *yi;
  double *yi_1d;
  double *zd;
  double *zdm;
  double *zi;

  mx = m;
  my = m;

  cout << "\n";
  cout << "LAGRANGE_INTERP_2D_TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_2D problem #" << prob << "\n";
  cout << "  Using polynomial interpolant of product degree " << mx << " x " << my << "\n";

  nd = ( mx + 1 ) * ( my + 1 );
  cout << "  Number of data points = " << nd << "\n";

  xd_1d = r8vec_chebyspace_new ( mx + 1, 0.0, 1.0 );
  yd_1d = r8vec_chebyspace_new ( my + 1, 0.0, 1.0 );

  xd = new double[(mx+1)*(my+1)];
  yd = new double[(mx+1)*(my+1)];
  zd = new double[(mx+1)*(my+1)];

  ij = 0;
  for ( j = 0; j < my + 1; j++ )
  {
    for ( i = 0; i < mx + 1; i++ )
    {
      xd[ij] = xd_1d[i];
      yd[ij] = yd_1d[j];
      ij = ij + 1;
    }
  }

  f00_f0 ( prob, nd, xd, yd, zd );

  if ( nd <= 20 )
  {
    r8vec3_print ( nd, xd, yd, zd, "  X, Y, Z data:" );
  }
//
//  #1:  Does interpolant match function at data points?
//
  ni = nd;

  xi = new double[ni];
  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    xi[i] = xd[i];
    yi[i] = yd[i];
  }

  zi = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );

  if ( ni <= 20 )
  {
    r8vec3_print ( ni, xi, yi, zi, "  X, Y, Z interpolation:" );
  }

  int_error = r8vec_norm_affine ( ni, zi, zd ) / ( double ) ( ni );

  cout << "\n";
  cout << "  RMS data interpolation error = " << int_error << "\n";

  delete [] xi;
  delete [] yi;
  delete [] zi;
//
//  #2:  Does interpolant approximate data at midpoints?
//
  if ( 1 < nd )
  {
    xi_1d = new double[mx];
    yi_1d = new double[my];

    for ( i = 0; i < mx; i++ )
    {
      xi_1d[i] = 0.5 * ( xd_1d[i] + xd_1d[i+1] );
    }
    for ( i = 0; i < my; i++ )
    {
      yi_1d[i] = 0.5 * ( yd_1d[i] + yd_1d[i+1] );
    }

    ni = mx * my;

    xi = new double[ni];
    yi = new double[ni];
    zdm = new double[ni];
    
    ij = 0;
    for ( j = 0; j < my; j++ )
    {
      for ( i = 0; i < mx; i++ )
      {
        xi[ij] = xi_1d[i];
        yi[ij] = yi_1d[j];
        ij = ij + 1;
      }
    }

    f00_f0 ( prob, ni, xi, yi, zdm );

    zi = lagrange_interp_2d ( mx, my, xd_1d, yd_1d, zd, ni, xi, yi );

    app_error = r8vec_norm_affine ( ni, zi, zdm ) / ( double ) ( ni );

    cout << "\n";
    cout << "  RMS data approximation error = " << app_error << "\n";

    delete [] xi;
    delete [] xi_1d;
    delete [] yi;
    delete [] yi_1d;
    delete [] zdm;
    delete [] zi;
  }

  delete [] xd;
  delete [] xd_1d;
  delete [] yd;
  delete [] yd_1d;
  delete [] zd;

  return;
}
