# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "laguerre_polynomial.hpp"

//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void imtqlx ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentialy) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cout << "\n";
        cout << "IMTQLX - Fatal error%\n";
        cout << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************80

double *l_exponential_product ( int p, double b )

//****************************************************************************80
//
//  Purpose:
//
//    L_EXPONENTIAL_PRODUCT: exponential product table for L(n,x).
//
//  Discussion:
//
//    Let L(n,x) represent the Laguerre polynomial of degree n.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of exp(B*X) with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( 0 <= X < +oo ) exp(b*x) * L(i,x) * L(j,x) * exp (-x) dx
//
//    Because of the exponential factor, the quadrature will not be exact.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, double B, the coefficient of X in the exponential factor.
//
//    Output, double L_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//    TABLE(I,J) represents the weighted integral of exp(B*X) * L(i,x) * L(j,x).
//
{
  int i;
  int j;
  int k;
  double *l_table;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = ( 3 * p + 4 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  l_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    l_table = l_polynomial ( 1, p, x_table+k );
    for ( j = 0; j <= p; j++ )
    {
      for ( i = 0; i <= p; i++ )
      {
        table[i+j*(p+1)] = table[i+j*(p+1)] 
          + w_table[k] * exp ( b * x ) * l_table[i] * l_table[j];
      }
    }
    delete [] l_table;
  }
  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

double l_integral ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    L_INTEGRAL evaluates a monomial integral associated with L(n,x).
//
//  Discussion:
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the exponent.
//    0 <= N.
//
//    Output, double L_INTEGRAL, the value of the integral.
//
{
  double value;

  value = r8_factorial ( n );

  return value;
}
//****************************************************************************80

double *l_polynomial ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    L_POLYNOMIAL evaluates the Laguerre polynomials L(n,x).
//
//  Differential equation:
//
//    X * Y'' + (1-X) * Y' + N * Y = 0
//
//  First terms:
//
//      1
//     -X    +  1
//   (  X^2 -  4 X     +  2 ) / 2
//   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
//   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
//   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
//   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
//   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
//      + 52920 X^2 - 35280 X + 5040 ) / 5040
//
//  Recursion:
//
//    L(0,X) = 1,
//    L(1,X) = 1-X,
//    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
//
//  Orthogonality:
//
//    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
//    = 0 if N /= M
//    = 1 if N == M
//
//  Special values:
//
//    L(N,0) = 1.
//
//  Relations:
//
//    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int M, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double X[M], the evaluation points.
//
//    Output, double L_POLYNOMIAL[M*(N+1)], the function values.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 1.0 - x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = ( ( ( double ) ( 2 * j - 1 ) - x[i] ) * v[i+(j-1)*m] 
                   + ( double ) (   - j + 1 )          * v[i+(j-2)*m] ) 
                   / ( double ) (     j     );
    }
  }

  return v;
}
//****************************************************************************80

double *l_polynomial_coefficients ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    L_POLYNOMIAL_COEFFICIENTS: coeffs for Laguerre polynomial L(n,x).
//
//  First terms:
//
//    0: 1
//    1: 1  -1
//    2: 1  -2  1/2
//    3: 1  -3  3/2  1/6
//    4: 1  -4  4   -2/3  1/24
//    5: 1  -5  5   -5/3  5/24  -1/120
//
//  Recursion:
//
//    L(0,X) = ( 1,  0, 0, ..., 0 )
//    L(1,X) = ( 1, -1, 0, ..., 0 )
//    L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Output, double L_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
//    of the Laguerre polynomials of degree 0 through N. 
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    c[i+0*(n+1)] = 1.0;
  }

  if ( n == 0 )
  {
    return c;
  }

  c[1+1*(n+1)] = -1.0;

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      c[i+j*(n+1)] = ( 
          ( double ) ( 2 * i - 1 ) * c[i-1+j*(n+1)] 
        + ( double ) (   - i + 1 ) * c[i-2+j*(n+1)] 
        -                            c[i-1+(j-1)*(n+1)] ) 
        / ( double )       i;
    }
  }
  return c;
}
//****************************************************************************80

void l_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    L_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial L(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,x]
//
//  Differential equation:
//
//    X * Y'' + (1-X) * Y' + N * Y = 0;
//
//  First terms:
//
//      1
//     -X    +  1
//   (  X^2 -  4 X     +  2 ) / 2
//   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
//   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
//   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
//   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
//   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
//      + 52920 X^2 - 35280 X + 5040 ) / 5040
//
//  Recursion:
//
//    L(0,X) = 1,
//    L(1,X) = 1-X,
//    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
//
//  Orthogonality:
//
//    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
//    = 0 if N /= M
//    = 1 if N == M
//
//  Special values:
//
//    L(N,0) = 1.
//
//  Relations:
//
//    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 17

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.0000000000000000E+00,
     -0.5000000000000000E+00,
     -0.6666666666666667E+00,
     -0.6250000000000000E+00,
     -0.4666666666666667E+00,
     -0.2569444444444444E+00,
     -0.4047619047619048E-01,
      0.1539930555555556E+00,
      0.3097442680776014E+00,
      0.4189459325396825E+00,
      0.4801341790925124E+00,
      0.4962122235082305E+00,
     -0.4455729166666667E+00,
      0.8500000000000000E+00,
     -0.3166666666666667E+01,
      0.3433333333333333E+02  };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5 };

  static double x_vec[N_MAX] = {
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     0.5E+00,
     3.0E+00,
     5.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *l_polynomial_zeros ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    L_POLYNOMIAL_ZEROS: zeros of the Laguerre polynomial L(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order of the polynomial.
//
//    Output, double L_POLYNOMIAL_ZEROS[N], the zeros.
//
{
  double *bj;
  int i;
  double *w;
  double *x;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = 1.0;
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( i + 1 );
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 2 * i + 1 );
  }

  w = new double[n];
  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  delete [] bj;
  delete [] w;

  return x;
}
//****************************************************************************80

double *l_power_product ( int p, int e )

//****************************************************************************80
//
//  Purpose:
//
//    L_POWER_PRODUCT: power product table for L(n,x).
//
//  Discussion:
//
//    Let L(n,x) represent the Laguerre polynomial of degree n.  
//
//    For polynomial chaos applications, it is of interest to know the
//    value of the integrals of products of X^E with every possible pair
//    of basis functions.  That is, we'd like to form
//
//      Tij = Integral ( 0 <= X < +oo ) x^e * L(i,x) * L(j,x) * exp (-x) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int P, the maximum degree of the polyonomial 
//    factors.  0 <= P.
//
//    Input, int E, the exponent of X.
//    0 <= E.
//
//    Output, double L_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
//    TABLE(I,J) represents the weighted integral of X^E * L(i,x) * L(j,x).
//
{
  int i;
  int j;
  int k;
  double *l_table;
  int order;
  double *table;
  double *w_table;
  double x;
  double *x_table;

  table = new double[(p+1)*(p+1)];

  for ( j = 0; j <= p; j++ )
  {
    for ( i = 0; i <= p; i++ )
    {
      table[i+j*(p+1)] = 0.0;
    }
  }

  order = p + 1 + ( e + 1 ) / 2;

  x_table = new double[order];
  w_table = new double[order];

  l_quadrature_rule ( order, x_table, w_table );

  for ( k = 0; k < order; k++ )
  {
    x = x_table[k];
    l_table = l_polynomial ( 1, p, x_table+k );

    if ( e == 0 )
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] + w_table[k] * l_table[i] * l_table[j];
        }
      }
    }
    else
    {
      for ( j = 0; j <= p; j++ )
      {
        for ( i = 0; i <= p; i++ )
        {
          table[i+j*(p+1)] = table[i+j*(p+1)] 
            + w_table[k] * pow ( x, e ) * l_table[i] * l_table[j];
        }
      }
    }
    delete [] l_table;
  }
  delete [] w_table;
  delete [] x_table;

  return table;
}
//****************************************************************************80

void l_quadrature_rule ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    L_QUADRATURE_RULE: Gauss-Laguerre quadrature based on L(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = 1.0;
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( i + 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 2 * i + 1 );
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double *lf_function ( int mm, int n, double alpha, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
//
//  Recursion:
//
//    Lf(0,ALPHA,X) = 1
//    Lf(1,ALPHA,X) = 1+ALPHA-X
//
//    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X) 
//                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
//
//  Restrictions:
//
//    -1 < ALPHA
//
//  Special values:
//
//    Lf(N,0,X) = L(N,X).
//    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
//
//  Norm:
//
//    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
//    = Gamma ( N + ALPHA + 1 ) / N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, double ALPHA, the parameter.  -1.0 < ALPHA.
//
//    Input, double X[MM], the evaluation points.
//
//    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[mm*(n+1)];

  for ( i = 0; i < mm; i++ )
  {
    v[i+0*mm] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < mm; i++ )
  {
    v[i+1*mm] = alpha + 1.0 - x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( (   alpha + ( double ) ( 2 * i - 1 ) - x[i] ) * v[i+(j-1)*mm] 
                  + ( - alpha + ( double ) (   - i + 1 )        ) * v[i+(j-2)*mm] ) 
                    /           ( double ) (     i     );
    }
  }

  return v;
}
//****************************************************************************80

void lf_function_values ( int &n_data, int &n, double &a, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LF_FUNCTION_VALUES: some values of the Laguerre function Lf(n,alpha,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,a,x]
//
//    The functions satisfy the following differential equation:
//
//      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0;
//
//    Function values can be generated by the recursion:
//
//      Lf(0,ALPHA,X) = 1
//      Lf(1,ALPHA,X) = 1+ALPHA-X
//
//      Lf(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * Lf(N-1,ALPHA,X)
//                     - (N-1+ALPHA) * Lf(N-2,ALPHA,X) ) / N
//
//    The parameter ALPHA is required to be greater than -1.
//
//    For ALPHA = 0, the generalized Laguerre function Lf(N,ALPHA,X)
//    is equal to the Laguerre polynomial L(N,X).
//
//    For ALPHA integral, the generalized Laguerre function
//    Lf(N,ALPHA,X) equals the associated Laguerre polynomial Lm(N,ALPHA,X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &A, the parameter.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.00E+00,
     0.25E+00,
     0.50E+00,
     0.75E+00,
     1.50E+00,
     2.50E+00,
     5.00E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00 };

  static double fx_vec[N_MAX] = {
      0.3726399739583333E-01,
      0.3494791666666667E+00,
      0.8710042317708333E+00,
      0.1672395833333333E+01,
      0.6657625325520833E+01,
      0.2395726725260417E+02,
      0.2031344319661458E+03,
      0.1284193996800000E+02,
      0.5359924801587302E+01,
      0.9204589064126984E+00,
     -0.1341585114857143E+01,
     -0.2119726307555556E+01,
     -0.1959193658349206E+01,
      0.1000000000000000E+01,
      0.5450000000000000E+01,
      0.1720125000000000E+02,
      0.4110393750000000E+02,
      0.8239745859375000E+02,
      0.1460179186171875E+03,
      0.2359204608298828E+03 };

  static int n_vec[N_MAX] = {
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     8,
     8,
     8,
     8,
     8,
     8,
     0,
     1,
     2,
     3,
     4,
     5,
     6 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.20E+00,
     0.40E+00,
     0.60E+00,
     0.80E+00,
     1.00E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *lf_function_zeros ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    LF_FUNCTION_ZEROS returns the zeros of Lf(n,alpha,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double ALPHA, the exponent of the X factor.
//    ALPHA must be nonnegative.
//
//    Output, double LF_FUNCTION_ZEROS[N], the zeros.
//
{
  double *bj;
  int i;
  double i_r8;
  double *w;
  double *x;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = r8_gamma ( alpha + 1.0 );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i );
    bj[i] = ( i_r8 + 1.0 ) * ( i_r8 + 1.0 + alpha );
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i );
    x[i] = ( double ) ( 2.0 * i_r8 + 1.0 + alpha );
  }

  w = new double[n];
  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;
  delete [] w;

  return x;
}
//****************************************************************************80

double lf_integral ( int n, double alpha )

//****************************************************************************80
//
//  Purpose:
//
//    LF_INTEGRAL evaluates a monomial integral associated with Lf(n,alpha,x).
//
//  Discussion:
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) x^n * x^alpha * exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the exponent.
//    0 <= N.
//
//    Input, double ALPHA, the exponent of X in the weight function.
//
//    Output, double LF_INTEGRAL, the value of the integral.
//
{
  double arg;
  double value;

  arg = alpha + ( double ) ( n + 1 );

  value = r8_gamma ( arg );

  return value;
}
//****************************************************************************80

void lf_quadrature_rule ( int n, double alpha, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LF_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lf(n,alpha,x);
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double ALPHA, the exponent of the X factor.
//    ALPHA must be nonnegative.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double i_r8;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = r8_gamma ( alpha + 1.0 );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i );
    bj[i] = ( i_r8 + 1.0 ) * ( i_r8 + 1.0 + alpha );
  }

  for ( i = 0; i < n; i++ )
  {
    i_r8 = ( double ) ( i );
    x[i] = ( double ) ( 2.0 * i_r8 + 1.0 + alpha );
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double lm_integral ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    LM_INTEGRAL evaluates a monomial integral associated with Lm(n,m,x).
//
//  Discussion:
//
//    The integral:
//
//      integral ( 0 <= x < +oo ) x^n * x^m * exp ( -x ) dx
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the exponent.
//    0 <= N.
//
//    Input, int M, the parameter.
//    0 <= M.
//
//    Output, double LM_INTEGRAL, the value of the integral.
//
{
  double value;

  value = r8_factorial ( n + m );

  return value;
}
//****************************************************************************80

double *lm_polynomial ( int mm, int n, int m, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LM_POLYNOMIAL evaluates Laguerre polynomials Lm(n,m,x).
//
//  First terms:
//
//    M = 0
//
//    Lm(0,0,X) =   1
//    Lm(1,0,X) =  -X   +  1
//    Lm(2,0,X) =   X^2 -  4 X   +  2
//    Lm(3,0,X) =  -X^3 +  9 X^2 -  18 X   +    6
//    Lm(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +     24
//    Lm(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x   +  120
//    Lm(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
//
//    M = 1
//
//    Lm(0,1,X) =    0
//    Lm(1,1,X) =   -1,
//    Lm(2,1,X) =    2 X - 4,
//    Lm(3,1,X) =   -3 X^2 + 18 X - 18,
//    Lm(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
//
//    M = 2
//
//    Lm(0,2,X) =    0
//    Lm(1,2,X) =    0,
//    Lm(2,2,X) =    2,
//    Lm(3,2,X) =   -6 X + 18,
//    Lm(4,2,X) =   12 X^2 - 96 X + 144
//
//    M = 3
//
//    Lm(0,3,X) =    0
//    Lm(1,3,X) =    0,
//    Lm(2,3,X) =    0,
//    Lm(3,3,X) =   -6,
//    Lm(4,3,X) =   24 X - 96
//
//    M = 4
//
//    Lm(0,4,X) =    0
//    Lm(1,4,X) =    0
//    Lm(2,4,X) =    0
//    Lm(3,4,X) =    0
//    Lm(4,4,X) =   24
//
//  Recursion:
//
//    Lm(0,M,X)   = 1 
//    Lm(1,M,X)   = (M+1-X)
//
//    if 2 <= N:
//
//      Lm(N,M,X)   = ( (M+2*N-1-X) * Lm(N-1,M,X) 
//                   +   (1-M-N)    * Lm(N-2,M,X) ) / N
//
//  Special values:
//
//    For M = 0, the associated Laguerre polynomials Lm(N,M,X) are equal 
//    to the Laguerre polynomials L(N,X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int MM, the number of evaluation points.
//
//    Input, int N, the highest order polynomial to compute.
//
//    Input, int M, the parameter.  M must be nonnegative.
//
//    Input, double X[MM], the evaluation points.
//
//    Output, double LM_POLYNOMIAL[MM*(N+1)], the function values.
//
{
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[mm*(n+1)];

  for ( j = 0; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }

  for ( i = 0; i < mm; i++ )
  {
    v[i+0*mm] = 1.0;
  }

  if ( n == 0 )
  {
    return v;
  }

  for ( i = 0; i < mm; i++ )
  {
    v[i+1*mm] = ( double ) ( m + 1 ) - x[i];
  }

  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( ( double ) (   m + 2 * j - 1 ) - x[i] ) * v[i+(j-1)*mm] 
                    + ( double ) ( - m     - j + 1 )          * v[i+(j-2)*mm] ) 
                    / ( double ) (           j     );
    }
  }

  return v;
}
//****************************************************************************80

double *lm_polynomial_coefficients ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    LM_POLYNOMIAL_COEFFICIENTS: coefficients of Laguerre polynomial Lm(n,m,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to compute.
//    Note that polynomials 0 through N will be computed.
//
//    Input, int M, the parameter.
//
//    Output, double LM_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
//    of the Laguerre polynomials of degree 0 through N. 
//
{
  double *c;
  int i;
  int j;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      c[i+j*(n+1)] = 0.0;
    }
  }

  c[0+0*(n+1)] = 1.0;

  if ( n == 0 )
  {
    return c;
  }

  c[1+0*(n+1)] = ( double ) ( m + 1 );
  c[1+1*(n+1)] = -1.0;

  for ( i = 2; i <= n; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      c[i+j*(n+1)] = ( 
          ( double ) (   m + 2 * i - 1 ) * c[i-1+j*(n+1)] 
        + ( double ) ( - m     - i + 1 ) * c[i-2+j*(n+1)] ) 
        / ( double )             i;
    }
    for ( j = 1; j <= i; j++ )
    {
      c[i+j*(n+1)] = c[i+j*(n+1)] - c[i-1+(j-1)*(n+1)] / ( double ) i;
    }
  }
  return c;
}
//****************************************************************************80

void lm_polynomial_values ( int &n_data, int &n, int &m, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LM_POLYNOMIAL_VALUES: some values of the Laguerre polynomial Lm(n,m,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,m,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, int &M, the parameter.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1500000000000000E+01,
     0.1625000000000000E+01,
     0.1479166666666667E+01,
     0.1148437500000000E+01,
     0.4586666666666667E+00,
     0.2878666666666667E+01,
     0.8098666666666667E+01,
     0.1711866666666667E+02,
     0.1045328776041667E+02,
     0.1329019368489583E+02,
     0.5622453647189670E+02,
     0.7484729341779436E+02,
     0.3238912982762806E+03,
     0.4426100000097533E+03,
     0.1936876572288250E+04 };

  static int m_vec[N_MAX] = {
    0, 0, 0, 0,
    0, 1, 1, 1,
    1, 0, 1, 2,
    3, 2, 2, 3,
    3, 4, 4, 5 };

  static int n_vec[N_MAX] = {
    1,  2,  3,  4,
    5,  1,  2,  3,
    4,  3,  3,  3,
    3,  4,  5,  6,
    7,  8,  9, 10 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *lm_polynomial_zeros ( int n, int m )

//****************************************************************************80
//
//  Purpose:
//
//    LM_POLYNOMIAL_ZEROS returns the zeros for Lm(n,m,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int M, the parameter.
//    0 <= M.
//
//    Output, double X[N], the zeros.
//
{
  double *bj;
  int i;
  double *w;
  double *x;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = r8_factorial ( m );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( i + 1 ) * ( i + 1 + m );
  }

  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 2 * i + 1 + m );
  }

  w = new double[n];
  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  delete [] bj;

  return x;
}
//****************************************************************************80

void lm_quadrature_rule ( int n, int m, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LM_QUADRATURE_RULE: Gauss-Laguerre quadrature rule for Lm(n,m,x);
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, int M, the parameter.
//    0 <= M.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  double *bj;
  int i;
  double zemu;
//
//  Define the zero-th moment.
//
  zemu = r8_factorial ( m );
//
//  Define the Jacobi matrix.
//
  bj = new double[n];
  for ( i = 0; i < n; i++ )
  {
    bj[i] = ( double ) ( i + 1 ) * ( i + 1 + m );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( 2 * i + 1 + m );
  }

  w[0] = sqrt ( zemu );
  for ( i = 1; i < n; i++ )
  {
    w[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( n, x, bj, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = w[i] * w[i];
  }

  delete [] bj;

  return;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL computes the factorial of N.
//
//  Discussion:
//
//    factorial ( N ) = product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//
//    Output, double R8_FACTORIAL, the factorial of N.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates Gamma(X) for an R8.
//
//  Discussion:
//
//    The C MATH library includes a function GAMMA ( X ) which should be
//    invoked instead of this function.
//
//    This routine calculates the gamma function for a real argument X.
//
//    Computation is based on an algorithm outlined in reference 1.
//    The program uses rational functions that approximate the gamma
//    function to at least 20 significant decimal digits.  Coefficients
//    for the approximation over the interval (1,2) are unpublished.
//    Those for the approximation for 12 <= X are from reference 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody,
//    An Overview of Software Development for Special Functions,
//    in Numerical Analysis Dundee, 1975,
//    edited by GA Watson,
//    Lecture Notes in Mathematics 506,
//    Springer, 1976.
//
//    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
//    Charles Mesztenyi, John Rice, Henry Thatcher,
//    Christoph Witzgall,
//    Computer Approximations,
//    Wiley, 1968,
//    LC: QA297.C64.
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double R8_GAMMA, the value of the function.
//
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  bool parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = false;
  fact = 1.0;
  n = 0;
  y = x;
//
//  Argument is negative.
//
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = true;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Argument is positive.
//
  if ( y < eps )
  {
//
//  Argument < EPS.
//
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
//
//  0.0 < argument < 1.0.
//
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
//
//  1.0 < argument < 12.0.
//  Reduce argument if necessary.
//
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
//
//  Evaluate approximation for 1.0 < argument < 2.0.
//
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
//
//  Adjust result for case  0.0 < argument < 1.0.
//
    if ( y1 < y )
    {
      res = res / y1;
    }
//
//  Adjust result for case 2.0 < argument < 12.0.
//
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
//
//  Evaluate for 12.0 <= argument.
//
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
//
//  Final adjustments and return.
//
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec2_print ( int n, double a1[], double a2[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
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
