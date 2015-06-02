# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_2d.hpp"
# include "test_interp_2d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int prob, int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PWL_INTERP_2D_PRB.
//
//  Discussion:
//
//    PWL_INTERP_2D_TEST tests the PWL_INTERP_2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int n_test[5] = { 2, 3, 4, 5, 9 };
  int n_test_num = 5;
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "PWL_INTERP_2D_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the PWL_INTERP_2D library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  The test needs the TEST_INTERP_2D library.\n";

  prob_num = f00_num ( );
//
//  Numerical tests.
//
  for ( prob = 1; prob <= prob_num; prob++ )
  {
    for ( i = 0; i < n_test_num; i++ )
    {
      n = n_test[i];
      test01 ( prob, n );
    }
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PWL_INTERP_2D_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int prob, int n )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_2D_TEST01 tests PWL_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int N, the grid size in each dimension.
//
{
  double app_error;
  int i;
  int ij;
  double int_error;
  int j;
  int nd;
  int ni;
  int nxd;
  int nyd;
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

  nxd = n;
  nyd = n;

  cout << "\n";
  cout << "PWL_INTERP_2D_TEST01:\n";
  cout << "  Interpolate data from TEST_INTERP_2D problem # " << prob << "\n";
  cout << "  Using polynomial interpolant of product degree " << nxd << " x " << nyd << "\n";

  nd = nxd * nyd;
  cout << "  Number of data points = " << nd << "\n";

  xd_1d = r8vec_linspace_new ( nxd, 0.0, 1.0 );
  yd_1d = r8vec_linspace_new ( nyd, 0.0, 1.0 );

  xd = new double[nxd*nyd];
  yd = new double[nxd*nyd];
  zd = new double[nxd*nyd];

  ij = 0;
  for ( j = 0; j < nyd; j++ )
  {
    for ( i = 0; i < nxd; i++ )
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
  xi = r8vec_copy_new ( ni, xd );
  yi = r8vec_copy_new ( ni, yd );

  zi = pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi, yi );

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
    xi_1d = new double[nxd-1];
    yi_1d = new double[nyd-1];

    for ( i = 0; i < nxd - 1; i++ )
    {
      xi_1d[i] = 0.5 * ( xd_1d[i] + xd_1d[i+1] );
    }
    for ( i = 0; i < nyd - 1; i++ )
    {
      yi_1d[i] = 0.5 * ( yd_1d[i] + yd_1d[i+1] );
    }

    ni = ( nxd - 1 ) * ( nyd - 1 );

    xi = new double[ni];
    yi = new double[ni];
    zdm = new double[ni];

    ij = 0;
    for ( j = 0; j < nyd - 1; j++ )
    {
      for ( i = 0; i < nxd - 1; i++ )
      {
        xi[ij] = xi_1d[i];
        yi[ij] = yi_1d[j];
        ij = ij + 1;
      }
    }

    f00_f0 ( prob, ni, xi, yi, zdm );

    zi = pwl_interp_2d ( nxd, nyd, xd_1d, yd_1d, zd, ni, xi, yi );

    app_error = r8vec_norm_affine ( ni, zi, zdm ) / ( double ) ( ni );

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
