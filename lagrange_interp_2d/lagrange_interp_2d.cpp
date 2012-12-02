# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_2d.hpp"

//****************************************************************************80

double lagrange_basis_function_1d ( int mx, double xd[], int i, double xi ) 

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_BASIS_FUNCTION_1D evaluates a 1D Lagrange basis function.
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
//    Input, int MX, the degree of the basis function.
//
//    Input, double XD[MX+1], the interpolation nodes.
//
//    Input, int I, the index of the basis function.
//    0 <= I <= MX.
//
//    Input, double XI, the evaluation point.
//
//    Output, double LAGRANGE_BASIS_FUNCTION_1D, the value of the I-th Lagrange 1D 
//    basis function for the nodes XD, evaluated at XI.
//
{
  int j;
  double yi;

  yi = 1.0;

  if ( xi != xd[i] )
  {
    for ( j = 0; j < mx + 1; j++ )
    {
      if ( j != i )
      {
        yi = yi * ( xi - xd[j] ) / ( xd[i] - xd[j] );
      }
    }
  }

  return yi;
}
//****************************************************************************80

double *lagrange_interp_2d ( int mx, int my, double xd_1d[], double yd_1d[], 
  double zd[], int ni, double xi[], double yi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_2D evaluates the Lagrange interpolant for a product grid.
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
//    Input, int MX, MY, the polynomial degree in X and Y.
//
//    Input, double XD_1D[MX+1], YD_1D[MY+1], the 1D data locations.
//
//    Input, double ZD[(MX+1)*(MY+1)], the 2D array of data values.
//
//    Input, int NI, the number of 2D interpolation points.
//
//    Input, double XI[NI], YI[NI], the 2D interpolation points.
//
//    Output, double LAGRANGE_INTERP_2D[NI], the interpolated values.
//
{
  int i;
  int j;
  int k;
  int l;
  double lx;
  double ly;
  double *zi;

  zi = new double[ni];

  for ( k = 0; k < ni; k++ )
  {
    l = 0;
    zi[k] = 0.0;
    for ( j = 0; j < my + 1; j++ )
    {
      for ( i = 0; i < mx + 1; i++ )
      {
        lx = lagrange_basis_function_1d ( mx, xd_1d, i, xi[k] );
        ly = lagrange_basis_function_1d ( my, yd_1d, j, yi[k] );
        zi[k] = zi[k] + zd[l] * lx * ly;
        l = l + 1;
      }
    }
  }
  return zi;
}
