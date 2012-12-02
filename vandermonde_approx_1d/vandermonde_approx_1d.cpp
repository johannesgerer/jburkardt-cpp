# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "vandermonde_approx_1d.hpp"
# include "qr_solve.hpp"

//****************************************************************************80

double *vandermonde_approx_1d_coef ( int n, int m, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_1D_COEF computes a 1D polynomial approximant.
//
//  Discussion:
//
//    We assume the approximating function has the form
//
//      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
//
//    We have n data values (x(i),y(i)) which must be approximated:
//
//      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
//
//    This can be cast as an Nx(M+1) linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
//      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
//      [ .................. ] [ ... ] = [ ... ]
//      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
//
//    In the typical case, N is greater than M+1 (we have more data and equations
//    than degrees of freedom) and so a least squares solution is appropriate,
//    in which case the computed polynomial will be a least squares approximant
//    to the data.
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
//    Input, int N, the number of data points.
//
//    Input, int M, the degree of the polynomial.
//
//    Input, double X[N], Y[N], the data values.
//
//    Output, double VANDERMONDE_APPROX_1D_COEF[M+1], the coefficients of 
//    the approximating polynomial.  C(0) is the constant term, and C(M) 
//    multiplies X^M.
//
{
  double *a;
  double *c;

  a = vandermonde_approx_1d_matrix ( n, m, x );

  c = qr_solve ( n, m + 1, a, y );

  delete [] a;

  return c;
}
//****************************************************************************80

double *vandermonde_approx_1d_matrix ( int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_1D_MATRIX computes a Vandermonde 1D approximation matrix.
//
//  Discussion:
//
//    We assume the approximant has the form
//
//      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
//
//    We have n data values (x(i),y(i)) which must be approximated:
//
//      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
//
//    This can be cast as an Nx(M+1) linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
//      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
//      [ .................. ] [ ... ] = [ ... ]
//      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
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
//    Input, int N, the number of data points.
//
//    Input, int M, the degree of the polynomial.
//
//    Input, double X(N), the data values.
//
//    Output, double VANDERMONDE_APPROX_1D_MATRIX[N*(M+1)], the Vandermonde matrix for X.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*(m+1)];

  for ( i = 0; i < n; i++ )
  {
    a[i+0*n] = 1.0;
  }

  for ( j = 1; j <= m; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = a[i+(j-1)*n] * x[i];
    }
  }

  return a;
}
