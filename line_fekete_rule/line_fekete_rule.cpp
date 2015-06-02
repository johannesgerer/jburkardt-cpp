# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "line_fekete_rule.hpp"
# include "qr_solve.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *cheby_van1 ( int m, double a, double b, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_VAN1 returns the CHEBY_VAN1 matrix.
//
//  Discussion:
//
//    Normally, the Chebyshev polynomials are defined on -1 <= XI <= +1.
//    Here, we assume the Chebyshev polynomials have been defined on the
//    interval A <= X <= B, using the mapping
//      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
//    so that
//      ChebyAB(A,B;X) = Cheby(XI).
//
//    if ( I == 1 ) then
//      V(1,1:N) = 1;
//    elseif ( I == 2 ) then
//      V(2,1:N) = XI(1:N);
//    else
//      V(I,1:N) = 2.0 * XI(1:N) * V(I-1,1:N) - V(I-2,1:N);
//
//  Example:
//
//    N = 5, A = -1, B = +1, X = ( 1, 2, 3, 4, 5 )
//
//    1  1   1    1    1
//    1  2   3    4    5
//    1  7  17   31   49
//    1 26  99  244  485
//    1 97 577 1921 4801
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Nicholas Higham,
//    Stability analysis of algorithms for solving confluent
//    Vandermonde-like systems,
//    SIAM Journal on Matrix Analysis and Applications,
//    Volume 11, 1990, pages 23-41.
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, double A, B, the interval.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, double X[N], the vector that defines the matrix.
//
//    Output, double CHEBY_VAN1[M*N], the matrix.
//
{
  int i;
  int j;
  double *v;
  double xi;

  v = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    xi = ( - ( b - x[j] ) + ( x[j] - a ) ) / ( b - a );
    for ( i = 0; i < m; i++ )
    {
      if ( i == 0 )
      {
        v[i+j*m] = 1.0;
      }
      else if ( i == 1 )
      {
        v[i+j*m] = xi;
      }
      else
      {
        v[i+j*m] = 2.0 * xi * v[i-1+j*m] - v[i-2+j*m];
      }
    }
  }
  return v;
}
//****************************************************************************80

double *legendre_van ( int m, double a, double b, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_VAN returns the LEGENDRE_VAN matrix.
//
//  Discussion:
//
//    The LEGENDRE_VAN matrix is the Legendre Vandermonde-like matrix.
//
//    Normally, the Legendre polynomials are defined on -1 <= XI <= +1.
//    Here, we assume the Legendre polynomials have been defined on the
//    interval A <= X <= B, using the mapping
//      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
//    so that
//      Lab(A,B;X) = L(XI).
//
//    if ( I = 1 ) then
//      V(1,1:N) = 1
//    else if ( I = 2 ) then
//      V(2,1:N) = XI(1:N)
//    else
//      V(I,1:N) = ( (2*I-1) * XI(1:N) * V(I-1,1:N) - (I-1)*V(I-2,1:N) ) / I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//
//    Input, double A, B, the interval.
//
//    Input, int N, the number of columns of the matrix.
//
//    Input, double X[N], the vector that defines the matrix.
//
//    Output, double LEGENDRE_VAN[M*N], the matrix.
//
{
  int i;
  int j;
  double *v;
  double xi;

  v = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    xi = ( - ( b - x[j] ) + ( x[j] - a ) ) / ( b - a );
    for ( i = 0; i < m; i++ )
    {
      if ( i == 0 )
      {
        v[i+j*m] = 1.0;
      }
      else if ( i == 1 )
      {
        v[i+j*m] = xi;
      }
      else
      {
        v[i+j*m] = ( ( double ) ( 2 * i - 1 ) * xi * v[i-1+j*m] +
                     ( double ) (   - i + 1 ) *      v[i-2+j*m] )
                   / ( double ) (     i );
      }
    }
  }

  return v;
}
//****************************************************************************80

void line_fekete_chebyshev ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_FEKETE_CHEBYSHEV: approximate Fekete points in an interval [A,B].
//
//  Discussion:
//
//    We use the Chebyshev basis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Len Bos, Norm Levenberg,
//    On the calculation of approximate Fekete points: the univariate case,
//    Electronic Transactions on Numerical Analysis, 
//    Volume 30, pages 377-397, 2008.
//   
//  Parameters:
//
//    Input, int M, the number of basis polynomials.
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the number of sample points.
//    M <= N.
//
//    Input, double X(N), the coordinates of the sample points.
//
//    Output, int &NF, the number of Fekete points.
//    If the computation is successful, NF = M.
//
//    Output, double XF(NF), the coordinates of the Fekete points.
//
//    Output, double WF(NF), the weights of the Fekete points.
//
{
  int i;
  int j;
  double *mom;
  double r8_pi = 3.141592653589793;
  double *v;
  double *w;

  if ( n < m )
  {
    cerr << "\n";
    cerr << "LINE_FEKETE_CHEBYSHEV - Fatal error!\n";
    cerr << "  N < M.\n";
    exit ( 1 );
  }
//
//  Compute the Chebyshev-Vandermonde matrix.
//
  v = cheby_van1 ( m, a, b, n, x );
//
//  MOM(I) = Integral ( A <= x <= B ) Tab(A,B,I;x) dx
//
  mom = new double[m];

  mom[0] = r8_pi * ( b - a ) / 2.0;
  for ( i = 1; i < m; i++ )
  {
    mom[i] = 0.0;
  }
//
//  Solve the system for the weights W.
//
  w = qr_solve ( m, n, v, mom );
//
//  Extract the data associated with the nonzero weights.
//
  nf = 0;
  for ( j = 0; j < n; j++)
  {
    if ( w[j] != 0.0 )
    {
      if ( nf < m )
      {
        xf[nf] = x[j];
        wf[nf] = w[j];
        nf = nf + 1;
      }
    }
  }

  delete [] mom;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void line_fekete_legendre ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_FEKETE_LEGENDRE: approximate Fekete points in an interval [A,B].
//
//  Discussion:
//
//    We use the uniform weight and the Legendre basis:
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Len Bos, Norm Levenberg,
//    On the calculation of approximate Fekete points: the univariate case,
//    Electronic Transactions on Numerical Analysis, 
//    Volume 30, pages 377-397, 2008.
//    
//  Parameters:
//
//    Input, int M, the number of basis polynomials.
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the number of sample points.
//    M <= N.
//
//    Input, double X(N), the coordinates of the sample points.
//
//    Output, int &NF, the number of Fekete points.
//    If the computation is successful, NF = M.
//
//    Output, double XF(NF), the coordinates of the Fekete points.
//
//    Output, double WF(NF), the weights of the Fekete points.
//
{
  int i;
  int j;
  double *mom;
  double *v;
  double *w;

  if ( n < m )
  {
    cerr << "\n";
    cerr << "LINE_FEKETE_LEGENDRE - Fatal error!\n";
    cerr << "  N < M.\n";
    exit ( 1 );
  }
//
//  Compute the Legendre-Vandermonde matrix.
//
  v = legendre_van ( m, a, b, n, x );
//
//  MOM(i) = integral ( A <= X <= B ) Lab(A,B,I;X) dx
//
  mom = new double[m];
  mom[0] = b - a;
  for ( i = 1; i < m; i++ )
  {
    mom[i] = 0.0;
  }
//
//  Solve the system for the weights W.
//
  w = qr_solve ( m, n, v, mom );
//
//  Extract the data associated with the nonzero weights.
//
  nf = 0;
  for ( j = 0; j < n; j++)
  {
    if ( w[j] != 0.0 )
    {
      if ( nf < m )
      {
        xf[nf] = x[j];
        wf[nf] = w[j];
        nf = nf + 1;
      }
    }
  }

  delete [] mom;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void line_fekete_monomial ( int m, double a, double b, int n, double x[], 
  int &nf, double xf[], double wf[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_FEKETE_MONOMIAL: approximate Fekete points in an interval [A,B].
//
//  Discussion:
//
//    We use the uniform weight and the monomial basis:
//
//      P(j) = x^(j-1)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alvise Sommariva, Marco Vianello,
//    Computing approximate Fekete points by QR factorizations of Vandermonde 
//    matrices,
//    Computers and Mathematics with Applications,
//    Volume 57, 2009, pages 1324-1336.
//    
//  Parameters:
//
//    Input, int M, the number of basis polynomials.
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the number of sample points.
//    M <= N.
//
//    Input, double X(N), the coordinates of the sample points.
//
//    Output, int &NF, the number of Fekete points.
//    If the computation is successful, NF = M.
//
//    Output, double XF(NF), the coordinates of the Fekete points.
//
//    Output, double WF(NF), the weights of the Fekete points.
//
{
  int i;
  int j;
  double *mom;
  double *v;
  double *w;

  if ( n < m )
  {
    cerr << "\n";
    cerr << "LINE_FEKETE_MONOMIAL - Fatal error!\n";
    cerr << "  N < M.\n";
    exit ( 1 );
  }
//
//  Form the moments.
//
  mom = line_monomial_moments ( a, b, m );
//
//  Form the rectangular Vandermonde matrix V for the polynomial basis.
//
  v = new double[m*n];
  for ( j = 0; j < n; j++ )
  {
    v[0+j*m] = 1.0;
    for ( i = 1; i < m; i++ )
    {
      v[i+j*m] = v[i-1+j*m] * x[j];
    }
  }
//
//  Solve the system for the weights W.
//
  w = qr_solve ( m, n, v, mom );
//
//  Extract the data associated with the nonzero weights.
//
  nf = 0;
  for ( j = 0; j < n; j++)
  {
    if ( w[j] != 0.0 )
    {
      if ( nf < m )
      {
        xf[nf] = x[j];
        wf[nf] = w[j];
        nf = nf + 1;
      }
    }
  }

  delete [] mom;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

double *line_monomial_moments ( double a, double b, int m )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_MONOMIAL_MOMENTS computes monomial moments in [A,B].
//
//  Discussion:
//
//    We use the uniform weight and the shifted and scaled monomial basis:
//
//      P(a,b,i;x) = xi(a,b;x)^(i-1)
//       xi(a,b;x) = ( - ( b - x ) + ( x - a ) ) / ( b - a )
//
//    The i-th moment is
//
//      mom(i) = integral ( a <= x <= b ) P(a,b,i;x) dx
//             = integral ( a <= x <= b ) xi(a,b;x)^(i-1) dx
//             = 0.5 * ( b - a ) * integral ( -1 <= xi <= +1 ) xi^(i-1) dxi
//             = 0.5 * ( b - a ) * xi^i / i | ( -1 <= xi <= +1 )
//             = 0.5 * ( b - a ) * ( 1 - (-1)^i ) / i
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int M, the number of basis polynomials.
//
//    Output, double LINE_MONOMIAL_MOMENTS[M], the moments.
//
{
  int i;
  double *mom;

  mom = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    mom[i] = ( b - a ) * ( double ) ( ( i + 1 ) % 2 ) / ( double ) ( i + 1 );
  }

  return mom;
}

