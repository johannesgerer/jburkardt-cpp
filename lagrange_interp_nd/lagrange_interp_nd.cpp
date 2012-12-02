# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_nd.hpp"

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
//    30 September 2012
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
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] =  cos ( ( double ) ( n - i ) * pi 
                    / ( double ) ( n - 1     ) );
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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

int i4vec_product ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRODUCT multiplies the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_PRODUCT = 24
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector
//
//    Output, int I4VEC_PRODUCT, the product of the entries of A.
//
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
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

double *lagrange_interp_nd_grid ( int m, int n_1d[], double a[], double b[], 
  int nd )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
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
//    Input, int N_1D[M], the order of the 1D rule to be used
//    in each dimension.
//
//    Input, double A[M], B[M], the lower and upper limits.
//
//    Input, int ND, the number of points in the product grid.
//
//    Output, double LAGRANGE_INTERP_ND_GRID[M*ND], the points at which data 
//    is to be sampled.
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
    n = n_1d[i];
    x_1d = cc_compute_points ( n );
    for ( j = 0; j < n; j++ )
    {
      x_1d[j] = 0.5 * ( ( 1.0 - x_1d[j] ) * a[i] 
                      + ( 1.0 + x_1d[j] ) * b[i] );
    }
    r8vec_direct_product ( i, n, x_1d, m, nd, xd );
    delete [] x_1d;
  }

  return xd;
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
//    Input, double A[M], B[M], the lower and upper limits.
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
    delete [] x_1d;
  }

  return xd;
}
//****************************************************************************80

int lagrange_interp_nd_size ( int m, int n_1d[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
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
//    Input, int N_1D[M], the order of the 1D rule to be used
//    in each dimension.
//
//    Output, int LAGRANGE_INTERP_ND_SIZE, the number of points in the product grid.
//
{
  int nd;
//
//  Determine the number of data points.
//
  nd = i4vec_product ( m, n_1d );

  return nd;
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

double *lagrange_interp_nd_value ( int m, int n_1d[], double a[], double b[], 
  int nd, double zd[], int ni, double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
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
//    Input, int N_1D[M], the order of the 1D rule to be used
//    in each dimension.
//
//    Input, double A[M], B[M], the lower and upper limits.
//
//    Input, int ND, the number of points in the product grid.
//
//    Input, double ZD[ND], the function evaluated at the points XD.
//
//    Input, int NI, the number of points at which the 
//    interpolant is to be evaluated.
//
//    Input, double XI[M*NI], the points at which the interpolant is 
//    to be evaluated.
//
//    Output, double LAGRANGE_INTERP_ND_VALUE[NI], the interpolant evaluated 
//    at the points XI.
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
      n = n_1d[i];
      x_1d = cc_compute_points ( n );
      for ( k = 0; k < n; k++ )
      {
        x_1d[k] = 0.5 * ( ( 1.0 - x_1d[k] ) * a[i] 
                        + ( 1.0 + x_1d[k] ) * b[i] );
      }
      value = lagrange_base_1d ( n, x_1d, 1, xi+i+j*m );
      r8vec_direct_product2 ( i, n, value, m, nd, w );
      delete [] value;
      delete [] x_1d;
    }
    zi[j] = r8vec_dot_product ( nd, w, zd );
  }

  delete [] w;

  return zi;
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
//    Input, double A[M], B[M], the lower and upper limits.
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
      delete [] value;
      delete [] x_1d;
    }
    zi[j] = r8vec_dot_product ( nd, w, zd );
  }

  delete [] w;

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
//    30 September 2012
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

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2,
//      ...,
//      item JK of 1D rule K.
//
//    In particular,
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1:
//      Order = 4
//      X(1:4) = ( 1, 2, 3, 4 )
//
//    Rule 2:
//      Order = 3
//      X(1:3) = ( 10, 20, 30 )
//
//    Rule 3:
//      Order = 2
//      X(1:2) = ( 100, 200 )
//
//    Product Rule:
//      Order = 24
//      X(1:24) =
//        ( 1, 10, 100 )
//        ( 2, 10, 100 )
//        ( 3, 10, 100 )
//        ( 4, 10, 100 )
//        ( 1, 20, 100 )
//        ( 2, 20, 100 )
//        ( 3, 20, 100 )
//        ( 4, 20, 100 )
//        ( 1, 30, 100 )
//        ( 2, 30, 100 )
//        ( 3, 30, 100 )
//        ( 4, 30, 100 )
//        ( 1, 10, 200 )
//        ( 2, 10, 200 )
//        ( 3, 10, 200 )
//        ( 4, 10, 200 )
//        ( 1, 20, 200 )
//        ( 2, 20, 200 )
//        ( 3, 20, 200 )
//        ( 4, 20, 200 )
//        ( 1, 30, 200 )
//        ( 2, 30, 200 )
//        ( 3, 30, 200 )
//        ( 4, 30, 200 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values
//    for factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double X[FACTOR_NUM*POINT_NUM], the elements of the
//    direct product, which are built up gradually.
//
//  Local Parameters:
//
//    Local, int START, the first location of a block of values to set.
//
//    Local, int CONTIG, the number of consecutive values to set.
//
//    Local, int SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, int REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( j = 0; j < point_num; j++ )
    {
      for ( i = 0; i < factor_num; i++ )
      {
        x[i+j*factor_num] = 0.0;
      }
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( i = 0; i < factor_order; i++ )
  {
    start = 0 + i * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( j = start; j < start + contig; j++ )
      {
        x[factor_index+j*factor_num] = factor_value[i];
      }
      start = start + skip;
    }
  }
  contig = contig * factor_order;

  return;
}
//****************************************************************************80

void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2,
//      ...,
//      item JK of 1D rule K.
//
//    In particular,
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task for the weights W.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1:
//      Order = 4
//      W(1:4) = ( 2, 3, 5, 7 )
//
//    Rule 2:
//      Order = 3
//      W(1:3) = ( 11, 13, 17 )
//
//    Rule 3:
//      Order = 2
//      W(1:2) = ( 19, 23 )
//
//    Product Rule:
//      Order = 24
//      W(1:24) =
//        ( 2 * 11 * 19 )
//        ( 3 * 11 * 19 )
//        ( 4 * 11 * 19 )
//        ( 7 * 11 * 19 )
//        ( 2 * 13 * 19 )
//        ( 3 * 13 * 19 )
//        ( 5 * 13 * 19 )
//        ( 7 * 13 * 19 )
//        ( 2 * 17 * 19 )
//        ( 3 * 17 * 19 )
//        ( 5 * 17 * 19 )
//        ( 7 * 17 * 19 )
//        ( 2 * 11 * 23 )
//        ( 3 * 11 * 23 )
//        ( 5 * 11 * 23 )
//        ( 7 * 11 * 23 )
//        ( 2 * 13 * 23 )
//        ( 3 * 13 * 23 )
//        ( 5 * 13 * 23 )
//        ( 7 * 13 * 23 )
//        ( 2 * 17 * 23 )
//        ( 3 * 17 * 23 )
//        ( 5 * 17 * 23 )
//        ( 7 * 17 * 23 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
//    factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double W[POINT_NUM], the elements of the
//    direct product, which are built up gradually.
//
//  Local Parameters:
//
//    Local, integer START, the first location of a block of values to set.
//
//    Local, integer CONTIG, the number of consecutive values to set.
//
//    Local, integer SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, integer REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
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

double r8vec_norm_affine ( int n, double v0[], double v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.
//
//  Discussion:
//
//    The affine vector L2 norm is defined as:
//
//      R8VEC_NORM_AFFINE(V0,V1)
//        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double V0[N], the base vector.
//
//    Input, double V1[N], the vector.
//
//    Output, double R8VEC_NORM_AFFINE, the affine L2 norm.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

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
