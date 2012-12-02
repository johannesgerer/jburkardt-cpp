# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "vandermonde_approx_2d.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

//****************************************************************************80

int triangle_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NUM returns the N-th triangular number.
//
//  Definition:
//
//    The N-th triangular number T(N) is formed by the sum of the first
//    N ints:
//
//      T(N) = sum ( 1 <= I <= N ) I
//
//    By convention, T(0) = 0.
//
//  Formula:
//
//    T(N) = ( N * ( N + 1 ) ) / 2
//
//  First Values:
//
//     0
//     1
//     3
//     6
//    10
//    15
//    21
//    28
//    36
//    45
//    55
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the index of the desired number, 
//    which must be at least 0.
//
//    Output, int TRIANGLE_NUM, the N-th triangular number.
//
{
  int value;

  value = ( n * ( n + 1 ) ) / 2;

  return value;
}
//****************************************************************************80

double *vandermonde_approx_2d_coef ( int n, int m, double x[], double y[], 
  double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_2D_COEF computes a 2D polynomial approximant.
//
//  Discussion:
//
//    We assume the approximating function has the form of a polynomial
//    in X and Y of total degree M.
//
//      p(x,y) = c00 
//             + c10 * x                + c01 *  y
//             + c20 * x^2   + c11 * xy + c02 * y^2
//             + ...
//             + cm0 * x^(m) + ...      + c0m * y^m.
//
//    If we let T(K) = the K-th triangular number 
//            = sum ( 1 <= I <= K ) I
//    then the number of coefficients in the above polynomial is T(M+1).
//
//    We have n data locations (x(i),y(i)) and values z(i) to approximate:
//
//      p(x(i),y(i)) = z(i)
//
//    This can be cast as an NxT(M+1) linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
//      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
//      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
//      [ ...................... ] [ ... ] = [ ... ]
//      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
//
//    In the typical case, N is greater than T(M+1) (we have more data and 
//    equations than degrees of freedom) and so a least squares solution is 
//    appropriate, in which case the computed polynomial will be a least squares
//    approximant to the data.
//
//    The polynomial defined by the T(M+1) coefficients C could be evaluated 
//    at the Nx2-vector x by the command
//
//      pval = r8poly_value_2d ( m, c, n, x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data points.
//
//    Input, int M, the maximum degree of the polynomial.
//
//    Input, double X[N], Y[N] the data locations.
//
//    Input, double Z[N], the data values.
//
//    Output, double VANDERMONDE_APPROX_2D_COEF[T(M+1)], the 
//    coefficients of the approximating polynomial.  
//
{
  double *a;
  double *c;
  int tm;

  tm = triangle_num ( m + 1 );

  a = vandermonde_approx_2d_matrix ( n, m, tm, x, y );

  c = qr_solve ( n, tm, a, z );

  delete [] a;

  return c;
}
//****************************************************************************80

double *vandermonde_approx_2d_matrix ( int n, int m, int tm, double x[], 
  double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_APPROX_2D_MATRIX computes a Vandermonde 2D approximation matrix.
//
//  Discussion:
//
//    We assume the approximating function has the form of a polynomial
//    in X and Y of total degree M.
//
//      p(x,y) = c00 
//             + c10 * x                + c01 * y
//             + c20 * x^2   + c11 * xy + c02 * y^2
//             + ...
//             + cm0 * x^(m) + ...      + c0m * y^m.
//
//    If we let T(K) = the K-th triangular number 
//            = sum ( 1 <= I <= K ) I
//    then the number of coefficients in the above polynomial is T(M+1).
//
//    We have n data locations (x(i),y(i)) and values z(i) to approximate:
//
//      p(x(i),y(i)) = z(i)
//
//    This can be cast as an NxT(M+1) linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
//      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
//      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
//      [ ...................... ] [ ... ] = [ ... ]
//      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2012
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
//    Input, int TM, the M+1st triangular number.
//
//    Input, double X[N], Y[N], the data locations.
//
//    Output, double VANDERMONDE_APPROX_2D_MATRIX[N*TM], the Vandermonde matrix for X.
//
{
  double *a;
  int ex;
  int ey;
  int i;
  int j;
  int s;

  a = new double[n*tm];
  j = 0;

  for ( s = 0; s <= m; s++ )
  {
    for ( ex = s; 0 <= ex; ex-- )
    {
      ey = s - ex;
      for ( i = 0; i < n; i++ )
      {
        a[i+j*n] = pow ( x[i], ex ) * pow ( y[i], ey ); 
      }
      j = j + 1;
    }
  }
 
  return a;
}
