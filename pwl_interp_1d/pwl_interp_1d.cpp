# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_1d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *pwl_basis_1d ( int nd, double xd[], int k, int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the data points.
//
//    Input, int K, the index of the desired basis function,
//    1 <= K <= ND.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double PW_BASIS_1D[NI], the basis function at the 
//    interpolation points.
//
{
  double *bk;
  int i;
  double t;

  bk = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    bk[i] = 0.0;
  }

  if ( nd == 1 )
  {
    for ( i = 0; i < ni; i++ )
    {
      bk[i] = 1.0;
    }
    return bk;
  }

  for ( i = 0; i < ni; i++ )
  {
    if ( k == 0 && xi[i] <= xd[k] )
    {
      t = ( xi[i] - xd[k] ) / ( xd[k+1] - xd[k] );
      bk[i] = 1.0 - t;
    }
    else if ( k == nd - 1 && xd[k] <= xi[i] )
    {
      t = ( xi[i] - xd[k-1] ) / ( xd[k] - xd[k-1] );
      bk[i] = t;
    }
    else if ( xd[k-1] < xi[i] && xi[i] <= xd[k] )
    {
      t = ( xi[i] - xd[k-1] ) / ( xd[k] - xd[k-1] );
      bk[i] = t;
    }
    else if ( xd[k] <= xi[i] && xi[i] < xd[k+1] )
    {
      t = ( xi[i] - xd[k] ) / ( xd[k+1] - xd[k] );
      bk[i] = 1.0 - t;
    }
  }
  return bk;
}
//****************************************************************************80

double *pwl_interp_1d ( int nd, double xd[], double yd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_1D evaluates the piecewise linear interpolant.
//
//  Discussion:
//
//    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
//    linear function which interpolates the data (XD(I),YD(I)) for I = 1
//    to ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2012
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
//    Output, double PWL_INTERP_1D[NI], the interpolated values.
//
{
  int i;
  int k;
  double t;
  double *yi;

  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    yi[i] = 0.0;
  }

  if ( nd == 1 )
  {
    for ( i = 0; i < ni; i++ )
    {
      yi[i] = yd[0];
    }
    return yi;
  }

  for ( i = 0; i < ni; i++ )
  {
    if ( xi[i] <= xd[0] )
    {
      t = ( xi[i] - xd[0] ) / ( xd[1] - xd[0] );
      yi[i] = ( 1.0 - t ) * yd[0] + t * yd[1];
    }
    else if ( xd[nd-1] <= xi[i] )
    {
      t = ( xi[i] - xd[nd-2] ) / ( xd[nd-1] - xd[nd-2] );
      yi[i] = ( 1.0 - t ) * yd[nd-2] + t * yd[nd-1];
    }
    else
    {
      for ( k = 1; k < nd; k++ )
      {
        if ( xd[k-1] <= xi[i] && xi[i] <= xd[k] )
        {
          t = ( xi[i] - xd[k-1] ) / ( xd[k] - xd[k-1] );
          yi[i] = ( 1.0 - t ) * yd[k-1] + t * yd[k];
          break;
        }
      }
    }
  }
  return yi;
}
