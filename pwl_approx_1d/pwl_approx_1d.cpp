# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "pwl_approx_1d.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *pwl_approx_1d ( int nd, double xd[], double yd[], int nc, double xc[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_APPROX_1D determines the control values for a PWL approximant.
//
//  Discussion:
//
//    The piecewise linear approximant is defined by NC control pairs 
//    (XC(I),YC(I)) and approximates ND data pairs (XD(I),YD(I)).
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
//    Input, int ND, the number of data points.
//    ND must be at least 1.
//
//    Input, double XD[ND], the data points.
//
//    Input, double YD[ND], the data values.
//
//    Input, int NC, the number of control points.
//    NC must be at least 1.
//
//    Input, double XC[NC], the control points.  Set these with a
//    command like 
//      xc = r8vec_linspace_new ( nc, xmin, xmax );
//
//    Output, double PWL_APPROX_1D[NC], the control values.
//
{
  double *a;
  double *yc;
//
//  Define the NDxNC linear system that determines the control values.
//
  a = pwl_approx_1d_matrix ( nd, xd, yd, nc, xc );
//
//  Solve the system.
//
  yc = qr_solve ( nd, nc, a, yd );

  free ( a );

  return yc;
}
//****************************************************************************80

double *pwl_approx_1d_matrix ( int nd, double xd[], double yd[], int nc, 
  double xc[] )

//****************************************************************************80
//
//  Purpose:
//
//    PWL_APPROX_1D_MATRIX returns the matrix for the PWL approximant controls.
//
//  Discussion:
//
//    The value of the piecewise linear approximant, using control points XC
//    and control values YC, evaluated at the point XD, can be represented by
//
//      YD = A * YC
//
//    where A is a matrix whose values depend on XC and XD.
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
//    Input, int ND, the number of data points.
//    ND must be at least 1.
//
//    Input, double XD[ND], the data points.
//
//    Input, double YD[ND], the data values.
//
//    Input, int NC, the number of control points.
//    NC must be at least 1.
//
//    Input, double XC[NC], the control points.
//
//    Output, double PWL_APPROX_1D_MATRIX[ND*NC], the matrix.
//
{
  double *a;
  int i;
  int j;
  int k;
  double t;

  a = new double[nd*nc];

  for ( j = 0; j < nc; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      a[i+j*nd] = 0.0;
    }
  }

  for ( i = 0; i < nd; i++ )
  {
    k = nc - 2;
    for ( j = 1; j < nc - 1; j++ )
    {
      if ( xd[i] < xc[j] )
      {
        k = j - 1;
        break;
      }
    }
    t = ( xd[i] - xc[k] ) / ( xc[k+1] - xc[k] );
    a[i+k*nd]     = 1.0 - t;
    a[i+(k+1)*nd] =       t;
  }

  return a;
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
//    09 October 2012
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
