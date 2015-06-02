# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "sparse_count.hpp"

//****************************************************************************80

int cc_se_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    CC_SE_SIZE: Clenshaw Curtis Slow Exponential Growth.
//
//  Discussion:
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This calculation is much faster than a previous method.  It simply
//    computes the number of new points that are added at each level in the
//    1D rule, and then counts the new points at a given DIM_NUM dimensional
//    level vector as the product of the new points added in each dimension.
//
//    This approach will work for nested families, and may be extensible
//    to other families, and to mixed rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int CC_SE_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int l;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int o;
  int p;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;
  new_1d[1] = 2;

  p = 3;
  o = 3;

  for ( l = 2; l <= level_max; l++ )
  {
    p = 2 * l + 1;
    if ( o < p )
    {
      new_1d[l] = o - 1;
      o = 2 * o - 1;
    }
    else
    {
      new_1d[l] = 0;
    }
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

int cfn_e_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    CFN_E_SIZE: Closed Fully Nested Exponential Growth.
//
//  Discussion:
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This calculation is much faster than a previous method.  It simply
//    computes the number of new points that are added at each level in the
//    1D rule, and then counts the new points at a given DIM_NUM dimensional
//    level vector as the product of the new points added in each dimension.
//
//    This approach will work for nested families, and may be extensible
//    to other families, and to mixed rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int CFN_E_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int j;
  int l;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;
  new_1d[1] = 2;

  j = 1;
  for ( l = 2; l <= level_max; l++ )
  {
    j = j * 2;
    new_1d[l] = j;
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

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
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

int f2_se_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    F2_SE_SIZE: Fejer Type 2 Slow Exponential Growth.
//
//  Discussion:
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This calculation is much faster than a previous method.  It simply
//    computes the number of new points that are added at each level in the
//    1D rule, and then counts the new points at a given DIM_NUM dimensional
//    level vector as the product of the new points added in each dimension.
//
//    This approach will work for nested families, and may be extensible
//    to other families, and to mixed rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int F2_SE_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int l;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int o;
  int p;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;

  p = 1;
  o = 1;

  for ( l = 1; l <= level_max; l++ )
  {
    p = 2 * l + 1;
    if ( o < p )
    {
      new_1d[l] = o + 1;
      o = 2 * o + 1;
    }
    else
    {
      new_1d[l] = 0;
    }
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

int gp_se_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    GP_SE_SIZE: Gauss Patterson Slow Exponential Growth.
//
//  Discussion:
//
//    The Gauss-Patterson-Slow family assumes that, for the underlying 1D
//    rules, a precision of 2*L+1 is needed at level L.  Therefore, the
//    lowest possible order Gauss-Patterson rule is chosen that will achieve
//    that precision.  This retains a combination of the advantages of
//    nestedness and high accuracy.
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      0 <= LEVEL <= LEVEL_MAX.
//
//    This calculation is much faster than a previous method.  It simply
//    computes the number of new points that are added at each level in the
//    1D rule, and then counts the new points at a given DIM_NUM dimensional
//    level vector as the product of the new points added in each dimension.
//
//    This approach will work for nested families, and may be extensible
//    to other families, and to mixed rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int GP_SE_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int j;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int o;
  int *order_1d;
  int p;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Count the points in the 1D rule.
//
  order_1d = new int[level_max+1];
  order_1d[0] = 1;
  for ( level = 1; level <= level_max; level++ )
  {
    p = 5;
    o = 3;
    while ( p < 2 * level + 1 )
    {
      p = 2 * p + 1;
      o = 2 * o + 1;
    }
    order_1d[level] = o;
  }
//
//  Count the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;
  for ( level = 1; level <= level_max; level++ )
  {
    new_1d[level] = order_1d[level] - order_1d[level-1];
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;
  delete [] order_1d;

  return point_num;
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
//    02 June 2007
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
//    Input, int N, K, are the values of N and K.
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

int ofn_e_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    OFN_E_SIZE: Open Fully Nested Exponential Growth.
//
//  Discussion:
//
//    This calculation assumes that an exponential growth rule is being used,
//    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int OFN_E_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int j;
  int l;
  int level;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 1;
  for ( l = 1; l <= level_max; l++ )
  {
    new_1d[l] = 2 * new_1d[l-1];
  }
//
//  Count the number of points by counting the number of new points 
//  associated with each level vector.
//
  level_1d = new int[dim_num];

  point_num = 0;

  for ( level = 0; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * new_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

int onn_e_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    ONN_E_SIZE Open Non-Nested Exponential Growth.
//
//  Discussion:
//
//    This calculation assumes that an exponential growth rule is being used,
//    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int ONN_E_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int j;
  int l;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int point_num;
  int t;
  int temp;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the 1D order vector.
//
  order_1d = new int[level_max+1];

  temp = 2;
  for ( l = 0; l <= level_max; l++ )
  {
    order_1d[l] = temp - 1;
    temp = temp * 2;
  }

  level_1d = new int[dim_num];

  point_num = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * order_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

int onn_l_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    ONN_L_SIZE Open Non-Nested Linear Growth.
//
//  Discussion:
//
//    This calculation assumes that a linear growth rule is being used,
//    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int ONN_L_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int j;
  int l;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int point_num;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the 1D order vector.
//
  order_1d = new int[level_max+1];

  for ( l = 0; l <= level_max; l++ )
  {
    order_1d[l] = 2 * l + 1;
  }

  level_1d = new int[dim_num];

  point_num = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  for ( level = level_min; level <= level_max; level++ )
  {
    more = false;
    h = 0;
    t = 0;

    for ( ; ;)
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );

      v = 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        v = v * order_1d[level_1d[dim]];
      }

      point_num = point_num + v;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

int own_e_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    OWN_E_SIZE: Open Weakly Nested Exponential Growth.
//
//  Discussion:
//
//    This calculation assumes that an exponential growth rule is being used,
//    that is, that the 1D rules have orders 1, 3, 7, 15, 31, and so on.
//
//    This calculation assumes that the 1D family of quadrature rules 
//    contains only one repeated point, presumably the value 0.0.
//    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
//    Generalized Gauss-Hermite rules.
//
//    The routine then counts the number of unique abscissas that will
//    be generated for a sparse grid of given dimension and level.
//
//    The computation is complicated.  It starts by counting just those
//    abscissas which have no 0.0 in them.  This is relatively easy, since
//    it is like counting the points in a sparse grid that uses open 
//    non-nested rules, but for which the order of each rule is reduced by 1.
//
//    Then we have to count the abscissas with one 0.0, two 0.0's and so
//    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
//    so for a particular number K of zeroes we only need to count the case
//    where the first K entries are zero, and multiply by C(DIM_NUM,K).
//
//    To count the number of entries with K zeroes, (and assuming 0 < K),
//    then, we essentially count the number of abscissas in an open 
//    non-nested rule as before, but modifed so that the minimum level is 0,
//    rather than LEVEL_MAX - DIM_NUM + 1.
//
//    I will mention that this was a rather difficult computation to
//    figure out!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 January 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int OWN_E_SIZE, the number of points in the grid.
//
{
  int dim;
  int dim_num2;
  int h;
  int j;
  int l;
  int level;
  int level_min;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int point_num2;
  int t;
  int temp;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 0;
  temp = 4;
  for ( l = 1; l <= level_max; l++ )
  {
    new_1d[l] = temp - 2;
    temp = temp * 2;
  }
//
//  Count the nonzero points in the full dimensional table with the usual
//  LEVEL_MIN restriction.
//
//  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
//  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
//  multiplying by the appropriate combinatorial coefficient.
//
  point_num = 0;

  for ( dim_num2 = dim_num; 0 <= dim_num2; dim_num2-- )
  {
    if ( dim_num2 == dim_num )
    {
      level_min = i4_max ( 0, level_max - dim_num + 1 );
    }
    else
    {
      level_min = 0;
    }

    if ( dim_num2 == 0 )
    {
      point_num2 = 1;
    }
    else
    {
      level_1d = new int[dim_num2];

      point_num2 = 0;

      for ( level = level_min; level <= level_max; level++ )
      {
        more = false;
        h = 0;
        t = 0;

        for ( ; ; )
        {
          comp_next ( level, dim_num2, level_1d, &more, &h, &t );

          v = 1;
          for ( dim = 0; dim < dim_num2; dim++ )
          {
            v = v * new_1d[level_1d[dim]];
          }

          point_num2 = point_num2 + v;

          if ( !more )
          {
            break;
          }
        }
      }
      delete [] level_1d;
    }
    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2;
  }
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

int own_l2_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    OWN_L2_SIZE: Open Weakly Nested Linear 2 Growth.
//
//  Discussion:
//
//    This calculation assumes that a linear growth rule of size 2 is being used,
//    that is, that the 1D rules have orders 1, 3, 5, 7, 9, and so on.
//
//    This calculation assumes that the 1D family of quadrature rules 
//    contains only one repeated point, presumably the value 0.0.
//    This assumption holds for Gauss-Legendre, Gauss-Hermite and 
//    Generalized Gauss-Hermite rules.
//
//    The routine then counts the number of unique abscissas that will
//    be generated for a sparse grid of given dimension and level.
//
//    The computation is complicated.  It starts by counting just those
//    abscissas which have no 0.0 in them.  This is relatively easy, since
//    it is like counting the points in a sparse grid that uses open 
//    non-nested rules, but for which the order of each rule is reduced by 1.
//
//    Then we have to count the abscissas with one 0.0, two 0.0's and so
//    on to DIM_NUM zeros.  We are assuming this is an isotropic grid,
//    so for a particular number K of zeroes we only need to count the case
//    where the first K entries are zero, and multiply by C(DIM_NUM,K).
//
//    To count the number of entries with K zeroes, (and assuming 0 < K),
//    then, we essentially count the number of abscissas in an open 
//    non-nested rule as before, but modifed so that the minimum level is 0,
//    rather than LEVEL_MAX - DIM_NUM + 1.
//
//    I will mention that this was a rather difficult computation to
//    figure out!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int OWN_L_SIZE, the number of points in the grid.
//
{
  int dim;
  int dim_num2;
  int h;
  int j;
  int l;
  int level;
  int level_min;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int point_num2;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  new_1d[0] = 0;
  for ( l = 1; l <= level_max; l++ )
  {
    new_1d[l] = 2 * l;
  }
//
//  Count the nonzero points in the full dimensional table with the usual
//  LEVEL_MIN restriction.
//
//  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
//  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
//  multiplying by the appropriate combinatorial coefficient.
//
  point_num = 0;

  for ( dim_num2 = dim_num; 0 <= dim_num2; dim_num2-- )
  {
    if ( dim_num2 == dim_num )
    {
      level_min = i4_max ( 0, level_max - dim_num + 1 );
    }
    else
    {
      level_min = 0;
    }

    if ( dim_num2 == 0 )
    {
      point_num2 = 1;
    }
    else
    {
      level_1d = new int[dim_num2];

      point_num2 = 0;

      for ( level = level_min; level <= level_max; level++ )
      {
        more = false;
        h = 0;
        t = 0;

        for ( ; ; )
        {
          comp_next ( level, dim_num2, level_1d, &more, &h, &t );

          v = 1;
          for ( dim = 0; dim < dim_num2; dim++ )
          {
            v = v * new_1d[level_1d[dim]];
          }

          point_num2 = point_num2 + v;

          if ( !more )
          {
            break;
          }
        }
      }
      delete [] level_1d;
    }
    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2;
  }
  delete [] new_1d;

  return point_num;
}
//****************************************************************************80

int own_o_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    OWN_O_SIZE: Open Weakly Nested Odd Growth.
//
//  Discussion:
//
//    This calculation assumes that an odd growth rule is being used,
//    that is, that the 1D rules have orders 1, 3, 3, 5, 5, 7, 7, 9, 9, 
//    and so on.
//
//    This repetition of rules is permissible because the sparse grid rule
//    only requires that the 1D rule of level L have precision at least 2*L+1.
//    This is achievable for Gaussian rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int OWN_O_SIZE, the number of points in the grid.
//
{
  int dim;
  int dim_num2;
  int h;
  int j;
  int l;
  int level;
  int level_min;
  int *level_1d;
  bool more;
  int *new_1d;
  int point_num;
  int point_num2;
  int t;
  int v;
//
//  Special case.
//
  if ( level_max < 0 )
  {
    point_num = 0;
    return point_num;
  }

  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }

  if ( dim_num == 1 )
  {
    point_num = 1 + 2 * ( ( level_max + 1 ) / 2 );
    return point_num;
  }
//
//  Construct the vector that counts the new points in the 1D rule.
//
  new_1d = new int[level_max+1];

  for ( l = 0; l <= level_max; l++ )
  {
    new_1d[l] = 0;
  }

  for ( l = 1; l <= level_max; l = l + 2 )
  {
    new_1d[l] = l + 1;
  }
//
//  Count the nonzero points in the full dimensional table with the usual
//  LEVEL_MIN restriction.
//
//  Then count the points with 1, 2, 3, ... DIM_NUM zeroes, by counting
//  the nonzero points in a DIM_NUM2 table, with LEVEL_MIN set to 0, and
//  multiplying by the appropriate combinatorial coefficient.
//
  point_num = 0;

  for ( dim_num2 = dim_num; 0 <= dim_num2; dim_num2-- )
  {
    if ( dim_num2 == dim_num )
    {
      level_min = i4_max ( 0, level_max - dim_num + 1 );
    }
    else
    {
      level_min = 0;
    }

    if ( dim_num2 == 0 )
    {
      point_num2 = 1;
    }
    else
    {
      level_1d = new int[dim_num2];

      point_num2 = 0;

      for ( level = level_min; level <= level_max; level++ )
      {
        more = false;
        h = 0;
        t = 0;

        for ( ; ; )
        {
          comp_next ( level, dim_num2, level_1d, &more, &h, &t );

          v = 1;
          for ( dim = 0; dim < dim_num2; dim++ )
          {
            v = v * new_1d[level_1d[dim]];
          }
 
          point_num2 = point_num2 + v;

          if ( !more )
          {
            break;
          }
        }
      }
      delete [] level_1d;
    }
    point_num = point_num + i4_choose ( dim_num, dim_num2 ) * point_num2;
  }

  delete [] new_1d;

  return point_num;
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
