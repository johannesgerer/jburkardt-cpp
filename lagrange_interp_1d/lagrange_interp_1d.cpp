# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_1d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *lagrange_basis_1d ( int nd, double xd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_BASIS_1D evaluates the Lagrange basis polynomials.
//
//  Discussion:
//
//    Given ND distinct abscissas, XD(1:ND),
//    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
//    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
//    other abscissas.
//
//    A formal representation is:
//
//      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
//       ( T - T(J) ) / ( T(I) - T(J) )
//
//    This routine accepts a set of NI values at which all the Lagrange
//    basis polynomials should be evaluated.
//
//    Given data values YD at each of the abscissas, the value of the
//    Lagrange interpolating polynomial at each of the interpolation points
//    is then simple to compute by matrix multiplication:
//
//      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2012
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
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double LAGRANGE_BASIS[NI*ND], the values
//    of the Lagrange basis polynomials at the interpolation points.
//
{
  int i;
  int j;
  int k;
  double *lb;
//
//  Evaluate the polynomial.
//
  lb = new double[ni*nd];

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < ni; i++ )
    {
      lb[i+j*ni] = 1.0;
    }
  }

  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < nd; j++ )
    {
      if ( j != i )
      {
        for ( k = 0; k < ni; k++ )
        {
          lb[k+i*ni] = lb[k+i*ni] * ( xi[k] - xd[j] ) / ( xd[i] - xd[j] );
        }
      }
    }
  }

  return lb;
}
//****************************************************************************80

double *lagrange_value_1d ( int nd, double xd[], double yd[], int ni, 
  double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
//
//  Discussion:
//
//    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
//    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
//    to ND.
//
//    The Lagrange interpolant can be constructed from the Lagrange basis
//    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange 
//    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree 
//    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
//
//    Given data values YD at each of the abscissas, the value of the
//    Lagrange interpolant may be written as
//
//      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2012
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
//    Output, double LAGRANGE_VALUE_1D[NI], the interpolated values.
//
{
  double *lb;
  double *yi;

  lb = lagrange_basis_1d ( nd, xd, ni, xi );

  yi = r8mat_mv_new ( ni, nd, lb, yd );

  delete [] lb;

  return yi;
}
