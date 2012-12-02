# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "vandermonde_interp_1d.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *vandermonde_interp_1d_coef ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_INTERP_1D_COEF computes a 1D polynomial interpolant.
//
//  Discussion:
//
//    We assume the interpolant has the form
//
//      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
//
//    We have n data values (x(i),y(i)) which must be interpolated:
//
//      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
//
//    This can be cast as an NxN linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
//      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
//      [ ...................... ] [ ... ] = [ ... ]
//      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
//
//    and if the x values are distinct, the system is theoretically
//    invertible, so we can retrieve the coefficient vector c and
//    evaluate the interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double X[N], Y[N], the data values.
//
//    Output, double VANDERMONDE_INTERP_1D_COEF[N], the coefficients of the 
//    interpolating polynomial.
//
{
  double *a;
  double *c;

  a = vandermonde_interp_1d_matrix ( n, x );

  c = qr_solve ( n, n, a, y );

  delete [] a;

  return c;
}
//****************************************************************************80

double *vandermonde_interp_1d_matrix ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_INTERP_1D_MATRIX computes a Vandermonde 1D interpolation matrix.
//
//  Discussion:
//
//    We assume the interpolant has the form
//
//      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
//
//    We have n data values (x(i),y(i)) which must be interpolated:
//
//      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
//
//    This can be cast as an NxN linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
//      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
//      [ ...................... ] [ ... ] = [ ... ]
//      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
//
//    and if the x values are distinct, the matrix A is theoretically
//    invertible (though in fact, generally badly conditioned).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, double X[N], the data values.
//
//    Output, double VANDERMONDE_INTERP_1D_MATRIX[N*N], the Vandermonde matrix for X.
//
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a[i+0*n] = 1.0;
  }
  for ( j = 1; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = a[i+(j-1)*n] * x[i];
    }
  }
  return a;
}
