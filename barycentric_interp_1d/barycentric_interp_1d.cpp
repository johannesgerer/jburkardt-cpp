# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "barycentric_interp_1d.hpp"
# include "r8lib.hpp"

//****************************************************************************

double *lagcheby1_interp_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

//****************************************************************************
//
//  Purpose:
//
//    LAGCHEBY1_INTERP_1D evaluates the Lagrange Chebyshev 1 interpolant.
//
//  Discussion:
//
//    The weight vector WD computed below is only valid if the data points
//    XD are, as expected, the Chebyshev Type 1 points for [-1,+1], or a linearly 
//    mapped version for [A,B].  The XD values may be computed by:
//
//      xd = r8vec_cheby1space ( nd, a, b )
//
//    for instance.
//
//    Thanks to John Ferrier for pointing out that DENOM and NUMER needed
//    to be initialized, 16 September 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
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
//    Output, double LAGCHEBY1_INTERP_1D[NI], the interpolated values.
//
{
  double *denom;
  int *exact;
  int i;
  int j;
  double *numer;
  double pi = 3.141592653589793;
  double t;
  double theta;
  double wd;
  double *yi;

  denom = new double[ni];
  exact = new int[ni];
  numer = new double[ni];
  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    exact[i] = -1;
    denom[i] = 0.0;
    numer[i] = 0.0;
  }

  for ( j = 0; j < nd; j++ )
  {
    theta = ( double ) ( 2 * j - 1 ) * pi / ( double ) ( 2 * nd );
    wd = r8_mop ( j + 1 ) * sin ( theta );

    for ( i = 0; i < ni; i++ )
    {
      if ( xi[i] == xd[j] )
      {
        exact[i] = j;
        numer[i] = yd[j];
        denom[i] = 1.0;
      }

      if ( exact[i] == -1 )
      {
        t = wd / ( xi[i] - xd[j] );
        numer[i] = numer[i] + t * yd[j];
        denom[i] = denom[i] + t;
      }
    }
  }

  for ( i = 0; i < ni; i++ )
  {
    yi[i] = numer[i] / denom[i];
  }

  return yi;
}
//****************************************************************************

double *lagcheby2_interp_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

//****************************************************************************
//
//  Purpose:
//
//    LAGCHEBY2_INTERP_1D evaluates the Lagrange Chebyshev 2 interpolant.
//
//  Discussion:
//
//    The weight vector WD computed below is only valid if the data points
//    XD are, as expected, the Chebyshev Type 2 points for [-1,+1], or a linearly 
//    mapped version for [A,B].  The XD values may be computed by:
//
//      xd = r8vec_cheby2space ( nd, a, b )
//
//    for instance.
//
//    Thanks to John Ferrier for pointing out that DENOM and NUMER needed
//    to be initialized, 16 September 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
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
//    Output, double LAGCHEBY2_INTERP_1D[NI], the interpolated values.
//
{
  double *denom;
  int *exact;
  int i;
  int j;
  double *numer;
  double t;
  double wd;
  double *yi;

  denom = new double[ni];
  exact = new int[ni];
  numer = new double[ni];
  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    exact[i] = -1;
    denom[i] = 0.0;
    numer[i] = 0.0;
  }

  for ( j = 0; j < nd; j++ )
  {
    wd = r8_mop ( j );
    if ( j == 0 || j == nd - 1 )
    {
      wd = 0.5 * wd;
    }

    for ( i = 0; i < ni; i++ )
    {
      if ( xi[i] == xd[j] )
      {
        exact[i] = j;
        numer[i] = yd[j];
        denom[i] = 1.0;
      }

      if ( exact[i] == -1 )
      {
        t = wd / ( xi[i] - xd[j] );
        numer[i] = numer[i] + t * yd[j];
        denom[i] = denom[i] + t;
      }
    }
  }

  for ( i = 0; i < ni; i++ )
  {
    yi[i] = numer[i] / denom[i];
  }

  return yi;
}
//****************************************************************************

double *lageven_interp_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

//****************************************************************************
//
//  Purpose:
//
//    LAGEVEN_VALUE_1D evaluates the Lagrange evenly-spaced interpolant.
//
//  Discussion:
//
//    The weight vector WD computed below is only valid if the data points
//    XD are, as expected, evenly spaced in an interval [A,B] with
//    spacing (B-A)/N.  The XD values might be computed by:
//
//      xd[i] = ( ( 2 * nd - 2 * i + 1 ) * a 
//              + (          2 * i - 1 ) * b ) 
//              / ( 2 * nd             )
//
//    for instance.
//
//    Thanks to John Ferrier for pointing out that DENOM and NUMER needed
//    to be initialized, 16 September 2013.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
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
//    Output, double LAGEVEN_INTERP_1D[NI], the interpolated values.
//
{
  double *denom;
  int *exact;
  int i;
  int j;
  double *numer;
  double t;
  double wd;
  double *yi;

  denom = new double[ni];
  exact = new int[ni];
  numer = new double[ni];
  yi = new double[ni];
  
  for ( i = 0; i < ni; i++ )
  {
    exact[i] = -1;
    denom[i] = 0.0;
    numer[i] = 0.0;
  }

  for ( j = 0; j < nd; j++ )
  {
    wd = r8_mop ( j ) * r8_choose ( nd, j );

    for ( i = 0; i < ni; i++ )
    {
      if ( xi[i] == xd[j] )
      {
        exact[i] = j;
        numer[i] = yd[j];
        denom[i] = 1.0;
      }

      if ( exact[i] == -1 )
      {
        t = wd / ( xi[i] - xd[j] );
        numer[i] = numer[i] + t * yd[j];
        denom[i] = denom[i] + t;
      }
    }
  }

  for ( i = 0; i < ni; i++ )
  {
    yi[i] = numer[i] / denom[i];
  }

  return yi;
}
