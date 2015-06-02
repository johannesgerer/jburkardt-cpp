# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "bernstein_polynomial.hpp"

//****************************************************************************80

double *bernstein_matrix ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_MATRIX returns the Bernstein matrix.
//
//  Discussion:
//
//    The Bernstein matrix of order N is an NxN matrix A which can be used to
//    transform a vector of power basis coefficients C representing a polynomial 
//    P(X) to a corresponding Bernstein basis coefficient vector B:
//
//      B = A * C
//
//    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
//    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
//    For N = 5, the matrix has the form:
//
//      1 -4   6  -4  1
//      0  4 -12  12 -4
//      0  0   6 -12  6
//      0  0   0   4 -4
//      0  0   0   0  1
//
//    and the numbers in each column represent the coefficients in the power
//    series expansion of a Bernstein polynomial, so that 
//
//      B(5,4) = - 4 x^4 + 12 x^3 - 12 x^2 + 4 x
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_MATRIX[N*N], the Bernstein matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];
 
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*n] = r8_mop ( j - i ) * r8_choose ( n - 1 - i, j - i ) 
        * r8_choose ( n - 1, i );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *bernstein_matrix_inverse ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
//
//  Discussion:
//
//    The inverse Bernstein matrix of order N is an NxN matrix A which can 
//    be used to transform a vector of Bernstein basis coefficients B
//    representing a polynomial P(X) to a corresponding power basis 
//    coefficient vector C:
//
//      C = A * B
//
//    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
//    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
//    For N = 5, the matrix has the form:
//
//      1   1    1    1   1
//      0  1/4  1/2  3/4  1
//      0   0   1/6  1/2  1
//      0   0    0   1/4  1
//      0   0    0    0   1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Output, double BERNSTEIN_MATRIX_INVERSE[N*N], the inverse Bernstein matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i <= j; i++ )
    {
      a[i+j*n] = r8_choose ( j, i ) / r8_choose ( n - 1, i );
    }
    for ( i = j + 1; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  return a;
}
//****************************************************************************80

double *bernstein_poly_01 ( int n, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01 evaluates the Bernstein polynomials based in [0,1].
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula is:
//
//      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials 
//    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
//    each of degree N, which form a basis for polynomials on [0,1].
//
//    Input, double X, the evaluation point.
//
//    Output, double BERNSTEIN_POLY[N+1], the values of the N+1 
//    Bernstein polynomials at X.
//
{
  double *bern;
  int i;
  int j;

  bern = new double[n+1];

  if ( n == 0 )
  {
    bern[0] = 1.0;
  }
  else if ( 0 < n )
  {
    bern[0] = 1.0 - x;
    bern[1] = x;
 
    for ( i = 2; i <= n; i++ )
    {
      bern[i] = x * bern[i-1];
      for ( j = i - 1; 1 <= j; j-- )
      {
        bern[j] =         x   * bern[j-1] 
                + ( 1.0 - x ) * bern[j];
      }
      bern[0] = ( 1.0 - x ) * bern[0];
    }
  }
  return bern;
}
//****************************************************************************80

void bernstein_poly_01_values ( int *n_data, int *n, int *k, double *x, 
  double *b )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_01_VALUES returns some values of the Bernstein polynomials.
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula for the Bernstein polynomials is
//
//      B(N,I)(X) = [N!/(I!(N-I)!)] * (1-X)^(N-I) * X^I
//
//    In Mathematica, the function can be evaluated by:
//
//      Binomial[n,i] * (1-x)^(n-i) * x^i
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2  * X
//    B(3,2)(X) = 3 * (1-X)    * X^2
//    B(3,3)(X) =                X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3  * X
//    B(4,2)(X) = 6 * (1-X)^2  * X^2
//    B(4,3)(X) = 4 * (1-X)    * X^3
//    B(4,4)(X) =                X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
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
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the degree of the polynomial.
//
//    Output, int *K, the index of the polynomial.
//
//    Output, double *X, the argument of the polynomial.
//
//    Output, double *B, the value of the polynomial B(N,K)(X).
//
{
# define N_MAX 15

  static double b_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.7500000000000000E+00,
     0.2500000000000000E+00,
     0.5625000000000000E+00,
     0.3750000000000000E+00,
     0.6250000000000000E-01,
     0.4218750000000000E+00,
     0.4218750000000000E+00,
     0.1406250000000000E+00,
     0.1562500000000000E-01,
     0.3164062500000000E+00,
     0.4218750000000000E+00,
     0.2109375000000000E+00,
     0.4687500000000000E-01,
     0.3906250000000000E-02 };

  static int k_vec[N_MAX] = {
    0,
    0, 1,
    0, 1, 2,
    0, 1, 2, 3,
    0, 1, 2, 3, 4 };

  static int n_vec[N_MAX] = {
    0,
    1, 1,
    2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4, 4 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *k = 0;
    *x = 0.0;
    *b = 0.0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *k = k_vec[*n_data-1];
    *x = x_vec[*n_data-1];
    *b = b_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double *bernstein_poly_ab ( int n, double a, double b, double x )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_AB evaluates at X the Bernstein polynomials based in [A,B].
//
//  Discussion:
//
//    The formula is:
//
//      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
//
//  First values:
//
//    B(0,0)(X) =   1
//
//    B(1,0)(X) = (      B-X                ) / (B-A)
//    B(1,1)(X) = (                 X-A     ) / (B-A)
//
//    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
//    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
//    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
//
//    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
//    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
//    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
//    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
//
//    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
//    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
//    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
//    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
//    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomials 
//    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
//    each of degree N, which form a basis for polynomials on [A,B].
//
//    Input, double A, B, the endpoints of the interval on which the
//    polynomials are to be based.  A and B should not be equal.
//
//    Input, double X, the point at which the polynomials 
//    are to be evaluated.
//
//    Output, double BERNSTEIN_POLY_AB[N+1], the values of the N+1
//    Bernstein polynomials at X.
//
{
  double *bern;
  int i;
  int j;

  if ( b == a )
  {
    cerr << "\n";
    cerr << "BERNSTEIN_POLY_AB - Fatal error!\n";
    cerr << "  A = B = " << a << "\n";
    exit ( 1 );
  }

  bern = new double[n+1];

  if ( n == 0 )
  {
    bern[0] = 1.0;
   }
  else if ( 0 < n )
  {
    bern[0] = ( b - x ) / ( b - a );
    bern[1] = ( x - a ) / ( b - a );
 
    for ( i = 2; i <= n; i++ )
    {
      bern[i] = ( x - a ) * bern[i-1] / ( b - a );
      for ( j = i - 1; 1 <= j; j-- )
      {
        bern[j] = ( ( b - x     ) * bern[j]
                  + (     x - a ) * bern[j-1] ) 
                  / ( b     - a );
      }
      bern[0] = ( b - x ) * bern[0] / ( b - a );
    }
  }
  return bern;
}
//****************************************************************************80

double *bernstein_poly_ab_approx ( int n, double a, double b, double ydata[], 
  int nval, double xval[] )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_AB_APPROX: Bernstein approximant to F(X) on [A,B].
//
//  Formula:
//
//    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
//
//    where
//
//      X(I) = ( ( N - I ) * A + I * B ) / N
//      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
//
//  Discussion:
//
//    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
//    not an interpolant; in other words, its value is not guaranteed to equal
//    that of F at any particular point.  However, for a fixed interval
//    [A,B], if we let N increase, the Bernstein polynomial converges
//    uniformly to F everywhere in [A,B], provided only that F is continuous.
//    Even if F is not continuous, but is bounded, the polynomial converges
//    pointwise to F(X) at all points of continuity.  On the other hand,
//    the convergence is quite slow compared to other interpolation
//    and approximation schemes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters:
//
//    Input, int N, the degree of the Bernstein polynomial
//    to be used.  N must be at least 0.
//
//    Input, double A, B, the endpoints of the interval on which the
//    approximant is based.  A and B should not be equal.
//
//    Input, double YDATA[N+1], the data values at N+1 equally
//    spaced points in [A,B].  If N = 0, then the evaluation point should
//    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
//    ( (N-I)*A + I*B ) / N ).
//
//    Input, int NVAL, the number of points at which the
//    approximant is to be evaluated.
//
//    Input, double XVAL[NVAL], the point at which the Bernstein 
//    polynomial approximant is to be evaluated.  The entries of XVAL do not 
//    have to lie in the interval [A,B].
//
//    Output, double BPAB_APPROX[NVAL], the values of the Bernstein 
//    polynomial approximant for F, based in [A,B], evaluated at XVAL.
//
{
  double *bvec;
  int i;
  double *yval;

  yval = new double[nval];

  for ( i = 0; i < nval; i++ )
  {
//
//  Evaluate the Bernstein basis polynomials at XVAL.
//
    bvec = bernstein_poly_ab ( n, a, b, xval[i] );
//
//  Now compute the sum of YDATA(I) * BVEC(I).
//
    yval[i] = r8vec_dot_product ( n + 1, ydata, bvec );

    delete [] bvec;
  }

  return yval;
}
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

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  double value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }
  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double r8mat_is_identity ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
//
//  Discussion:
//
//    An R8MAT is a matrix of real ( kind = 8 ) values.
//
//    The routine returns the Frobenius norm of A - I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
//    of the difference matrix A - I, which would be exactly zero
//    if A were the identity matrix.
//
{
  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        t = a[i+j*n] - 1.0;
      }
      else
      {
        t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t * t;
    }
  }
  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
}
//****************************************************************************80

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

double r8mat_norm_fro ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
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
//    Input, double A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

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

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
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
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
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
