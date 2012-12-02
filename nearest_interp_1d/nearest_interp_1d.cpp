# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "nearest_interp_1d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *nearest_interp_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    NEAREST_INTERP_1D evaluates the nearest neighbor interpolant.
//
//  Discussion:
//
//    The nearest neighbor interpolant L(ND,XD,YD)(X) is the piecewise
//    constant function which interpolates the data (XD(I),YD(I)) for I = 1
//    to ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//    ND must be at least 1.
//
//    Input, double XD[ND], the data points.
//
//    Input, double YD[ND], the data values.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double NEAREST_INTERP_1D[NI], the interpolated values.
//
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double *yi;

  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    k = 0;
    d = r8_abs ( xi[i] - xd[k] );
    for ( j = 1; j < nd; j++ )
    {
      d2 = r8_abs ( xi[i] - xd[j] );
      if ( d2 < d )
      {
        k = j;
        d = d2;
      }
    }
    yi[i] = yd[k];
  }

  return yi;
}
