# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "chebyshev.hpp"

//****************************************************************************80

double *chebyshev_coefficients ( double a, double b, int n, 
  double f ( double x ) )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_COEFFICIENTS determines Chebyshev interpolation coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Second Edition,
//    Cambridge University Press, 1992,
//    ISBN: 0-521-43064-X,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, double A, B, the domain of definition.
//
//    Input, int N, the order of the interpolant.
//
//    Input, double F ( double X ), an external function.
//
//    Output, double CHEBYSHEV_COEFFICIENTS[N], the Chebyshev coefficients.
//
{
  double angle;
  double *c;
  double *fx;
  int i;
  int j;
  double pi = 3.141592653589793;
  double x;

  fx = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( 2 * i + 1 ) * pi / ( double ) ( 2 * n );
    x = cos ( angle );
    x = 0.5 * ( a + b ) + x * 0.5 * ( b - a );
    fx[i] = f ( x );
  }

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      angle = ( double ) ( i * ( 2 * j + 1 ) ) * pi / ( double ) ( 2 * n );
      c[i] = c[i] + fx[j] * cos ( angle );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    c[i] = 2.0 * c[i] / ( double ) ( n );
  }

  delete [] fx;
  
  return c;
}
//****************************************************************************80

double *chebyshev_interpolant ( double a, double b, int n, double c[], int m, 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_INTERPOLANT evaluates a Chebyshev interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Second Edition,
//    Cambridge University Press, 1992,
//    ISBN: 0-521-43064-X,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, double A, B, the domain of definition.
//
//    Input, int N, the order of the polynomial.
//
//    Input, double C[N], the Chebyshev coefficients.
//
//    Input, int M, the number of points.
//
//    Input, double X[M], the point at which the polynomial is
//    to be evaluated.
//
//    Output, double CHEBYSHEF_INTERPOLANT[M], the value of the Chebyshev
//    polynomial at X.
//
{
  double *cf;
  double di;
  double dip1;
  double dip2;
  int i;
  int j;
  double y;

  cf = new double[m];

  for ( j = 0; j < m; j++ )
  {
    dip1 = 0.0;
    di = 0.0;
    y = ( 2.0 * x[j] - a  - b ) / ( b - a );

    for ( i = n - 1; 1 <= i; i-- )
    {
      dip2 = dip1;
      dip1 = di;
      di = 2.0 * y * dip1 - dip2 + c[i];
    }
    cf[j] = y * di - dip1 + 0.5 * c[0];
  }

  return cf;
}
//****************************************************************************80

double *chebyshev_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV_ZEROS returns zeroes of the Chebyshev polynomial T(N)(X).
//
//  Discussion:
//
//    We produce the Chebyshev zeros in ascending order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double CHEBYSHEV_ZEROS[N], the zeroes of T(N)(X).
//
{
  double angle;
  int i;
  double pi = 3.141592653589793;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double) ( 2 * ( n - i ) - 1 ) * pi / ( double ) ( 2 * n );
    x[i] = cos ( angle );
  }

  return x;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
