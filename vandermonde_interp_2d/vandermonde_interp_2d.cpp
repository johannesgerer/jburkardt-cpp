# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "vandermonde_interp_2d.hpp"

//****************************************************************************80

int triangle_num ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NUM returns the N-th triangular number.
//
//  Discussion:
//
//    The N-th triangular number T(N) is formed by the sum of the first
//    N integers:
//
//      T(N) = sum ( 1 <= I <= N ) I
//
//    By convention, T(0) = 0.
//
//    T(N) can be computed quickly by the formula:
//
//      T(N) = ( N * ( N + 1 ) ) / 2
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
//    07 October 2012
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

double *vandermonde_interp_2d_matrix ( int n, int m, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    VANDERMONDE_INTERP_2D_MATRIX computes a Vandermonde 2D interpolation matrix.
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
//    and we assume that N = T(M+1).
//
//    This can be cast as an NxN linear system for the polynomial
//    coefficients:
//
//      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
//      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
//      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
//      [ ...................... ] [ ... ] = [ ... ]
//      [ 1 xn yn  xn^2 ... yn^m ] [ c0n ] = [  zn ]
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
//    Input, int N, the number of data points.  It is necessary 
//    that N = T(M+1), where T(K) is the K-th triangular number.
//
//    Input, int M, the degree of the polynomial.
//
//    Input, double X[N], Y[N], the data locations.
//
//    Output, double VANDERMONDE_INTERP_2D_MATRIX[N*N], the Vandermonde matrix for X.
//
{
  double *a;
  int ex;
  int ey;
  int i;
  int j;
  int s;
  int tmp1;

  tmp1 = triangle_num ( m + 1 );

  if ( n != tmp1 )
  {
    cerr << "\n";
    cerr << "VANDERMONDE_INTERP_2D_MATRIX - Fatal error!\n";
    cerr << "  For interpolation, we need N = T(M+1).\n";
    cerr << "  But we have N = " << n << "\n";
    cerr << "  M = " << m << "\n";
    cerr << "  and T(M+1) = " << tmp1 << "\n";
    exit ( 1 );
  }

  a = new double[n*n];
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
