# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <cmath>
# include <ctime>

using namespace std;

# include "sparse_interp_nd.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *cc_compute_points ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_COMPUTE_POINTS computes Clenshaw Curtis quadrature points.
//
//  Discussion:
//
//    Our convention is that the abscissas are numbered from left to right.
//
//    This rule is defined on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double CC_COMPUTE_POINTS[N], the abscissas.
//
{
  int i;
  double pi = 3.141592653589793;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CC_COMPUTE_POINTS - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] =  cos ( ( double ) ( n - 1 - i ) * pi 
                  / ( double ) ( n - 1 ) );
    }
    x[0] = -1.0;
    if ( ( n % 2 ) == 1 )
    {
      x[(n-1)/2] = 0.0;
    }
    x[n-1] = +1.0;
  }
  return x;
}
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool &more, int &h, int &t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool &MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int &H, &T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( more ) )
  {
    t = n;
    h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < t )
    {
      h = 0;
    }
    h = h + 1;
    t = a[h-1];
    a[h-1] = 0;
    a[0] = t - 1;
    a[h] = a[h] + 1;
  }

  more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

int i4_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
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
//    09 November 2007
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
//    Output, int I4_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

int i4_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MOP returns the I-th power of -1 as an I4 value.
//
//  Discussion:
//
//    An I4 is an int value.
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
//    Output, int I4_MOP, the I-th power of -1.
//
{
  int value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

double *lagrange_base_1d ( int nd, double xd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_BASE_1D evaluates the Lagrange basis polynomials.
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
//    13 September 2012
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

double *lagrange_interp_nd_grid2 ( int m, int ind[], double a[], double b[], 
  int nd )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int IND[M], the index or level of the 1D rule 
//    to be used in each dimension.
//
//    Input, double A[M], B[M], the upper and lower limits.
//
//    Input, int ND, the number of points in the product grid.
//
//    Output, double LAGRANGE_INTERP_ND_GRID2[M*ND], the points at which data 
//    was sampled.
//
{
  int i;
  int j;
  int n;
  double *x_1d;
  double *xd;
//
//  Compute the data points.
//
  xd = new double[m*nd];

  for ( j = 0; j < nd; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      xd[i+j*m] = 0.0;
    }
  }

  for ( i = 0; i < m; i++ )
  {
    n = order_from_level_135 ( ind[i] );
    x_1d = cc_compute_points ( n );
    for ( j = 0; j < n; j++ )
    {
      x_1d[j] = 0.5 * ( ( 1.0 - x_1d[j] ) * a[i] 
                      + ( 1.0 + x_1d[j] ) * b[i] );
    }
    r8vec_direct_product ( i, n, x_1d, m, nd, xd );
    free ( x_1d );
  }

  return xd;
}
//****************************************************************************80

int lagrange_interp_nd_size2 ( int m, int ind[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int IND[M], the index or level of the 1D rule 
//    to be used in each dimension.
//
//    Output, int LAGRANGE_INTERP_ND_SIZE2, the number of points in the product grid.
//
{
  int i;
  int n;
  int nd;
//
//  Determine the number of data points.
//
  nd = 1;
  for ( i = 0; i < m; i++ )
  {
    n = order_from_level_135 ( ind[i] );
    nd = nd * n;
  }

  return nd;
}
//****************************************************************************80

double *lagrange_interp_nd_value2 ( int m, int ind[], double a[], double b[], 
  int nd, double zd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int IND[M], the index or level of the 1D rule 
//    to be used in each dimension.
//
//    Input, double A[M], B[M], the upper and lower limits.
//
//    Input, int ND, the number of points in the product grid.
//
//    Input, double ZD[ND], the function evaluated at the points XD.
//
//    Input, int NI, the number of points at which the 
//    interpolant is to be evaluated.
//
//    Input, double XI[M*NI], the points at which the interpolant 
//    is to be evaluated.
//
//    Output, double ZI[NI], the interpolant evaluated at the 
//    points XI.
//
{
  int i;
  int j;
  int k;
  int n;
  double *value;
  double *w;
  double *x_1d;
  double *zi;

  w = new double[nd];
  zi = new double[ni];

  for ( j = 0; j < ni; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      w[i] = 1.0;
    }

    for ( i = 0; i < m; i++ )
    {
      n = order_from_level_135 ( ind[i] );
      x_1d = cc_compute_points ( n );
      for ( k = 0; k < n; k++ )
      {
        x_1d[k] = 0.5 * ( ( 1.0 - x_1d[k] ) * a[i] 
                        + ( 1.0 + x_1d[k] ) * b[i] );
      }
      value = lagrange_base_1d ( n, x_1d, 1, xi+i+j*m );
      r8vec_direct_product2 ( i, n, value, m, nd, w );
      free ( value );
      free ( x_1d );
    }
    zi[j] = r8vec_dot_product ( nd, w, zd );
  }

  free ( w );

  return zi;
}
//****************************************************************************80

int order_from_level_135 ( int l )

//****************************************************************************80
//
//  Purpose:
//
//    ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
//
//  Discussion:
//
//    Clenshaw Curtis rules, and some others, often use the following
//    scheme:
//
//    L: 0  1  2  3   4   5
//    N: 1  3  5  9  17  33 ... 2^L+1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, the level, which should be 0 or greater.
//
//    Output, int ORDER_FROM_LEVEL_135, the order.
//
{
  int n;

  if ( l < 0 )
  {
    cerr << "\n";
    cerr << "ORDER_FROM_LEVEL_135 - Fatal error!\n";
    cerr << "  Illegal input value of L!\n";
    exit ( 1 );
  }
  else if ( l == 0 )
  {
    n = 1;
  }
  else
  {
    n = i4_power ( 2, l ) + 1;
  }
  return n;
}
//****************************************************************************80

void smolyak_coefficients ( int l_max, int m, int c[], int w[] )

//****************************************************************************80
//
//  Purpose:
//
//    SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
//
//  Discussion:
//
//    The Smolyak sparse interpolant can be written as:
//
//      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
//        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
//
//    where:
//
//    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
//    * |L| is the sum of the entries of L;
//    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
//    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
//
//    Note that:
//
//    * W(|L|) will represent the number of distinct interpolants for which
//      the sublevel, or sum of the L vector entries, is |L|;
//
//    * the coefficients C and counts W will be zero for sublevels 
//      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
//
//    * it will be the case that W' * C = 1, essentially because the interpolant
//      to the identity function must be the identity function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L_MAX, the (maximum) level.
//    0 <= L_MA.
//
//    Input, int M, the spatial dimension.
//    1 <= M.
//
//    Output, int C[L_MAX+1], the coefficients for objects 
//    at sublevels 0 through L_MAX.
//
//    Output, int W[L_MAX+1], the number of objects at 
//    sublevels 0 through L_MAX.
//
{
  int l;
  int l_min;

  l_min = i4_max ( l_max - m + 1, 0 );

  for ( l = 0; l < l_min; l++ )
  {
    c[l] = 0;
  }
  for ( l = l_min; l <= l_max; l++ )
  {
    c[l] = i4_mop ( l_max - l ) * i4_choose ( m - 1, l_max - l );
  }

  for ( l = 0; l < l_min; l++ )
  {
    w[l] = 0;
  }
  for ( l = l_min; l <= l_max; l++ )
  {
    w[l] = i4_choose ( l + m - 1, m - 1 );
  }
  return;
}

