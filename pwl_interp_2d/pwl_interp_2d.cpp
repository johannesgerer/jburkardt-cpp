# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "pwl_interp_2d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *pwl_interp_2d ( int nxd, int nyd, double xd[], double yd[], double zd[], 
  int ni, double xi[], double yi[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
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
//    Input, int NXD, NYD, the number of X and Y data values.
//
//    Input, double XD[NXD], YD[NYD], the sorted X and Y data.
//
//    Input, double ZD[NXD*NYD}, the Z data.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], YI[NI], the coordinates of the
//    interpolation points.
//
//    Output, double PWL_INTERP_2D[NI], the value of the interpolant.
//
{
  double alpha;
  double beta;
  double det;
  double dxa;
  double dxb;
  double dxi;
  double dya;
  double dyb;
  double dyi;
  double gamma;
  int i;
  int j;
  int k;
  double *zi;

  zi = new double[ni];

  for ( k = 0; k < ni; k++ )
  {
    i = r8vec_bracket5 ( nxd, xd, xi[k] );
    if ( i == -1 )
    {
      zi[k] = r8_huge ( );
      continue;
    }

    j = r8vec_bracket5 ( nyd, yd, yi[k] );
    if ( j == -1 )
    {
      zi[k] = r8_huge ( );
      continue;
    }

    if ( yi[k] < yd[j+1] + ( yd[j] - yd[j+1] ) * ( xi[i] - xd[i] ) / ( xd[i+1] - xd[i] ) )
    {
      dxa = xd[i+1] - xd[i];
      dya = yd[j]   - yd[j];

      dxb = xd[i]   - xd[i];
      dyb = yd[j+1] - yd[j];

      dxi = xi[k]   - xd[i];
      dyi = yi[k]   - yd[j];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

      zi[k] = alpha * zd[i+1+j*nxd] + beta * zd[i+(j+1)*nxd] + gamma * zd[i+j*nxd];
    }
    else
    {
      dxa = xd[i]   - xd[i+1];
      dya = yd[j+1] - yd[j+1];

      dxb = xd[i+1] - xd[i+1];
      dyb = yd[j]   - yd[j+1];

      dxi = xi[k]   - xd[i+1];
      dyi = yi[k]   - yd[j+1];

      det = dxa * dyb - dya * dxb;

      alpha = ( dxi * dyb - dyi * dxb ) / det;
      beta =  ( dxa * dyi - dya * dxi ) / det;
      gamma = 1.0 - alpha - beta;

      zi[k] = alpha * zd[i+(j+1)*nxd] + beta * zd[i+1+j*nxd] + gamma * zd[i+1+(j+1)*nxd];
    }
  }

  return zi;
}
