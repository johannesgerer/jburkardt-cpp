# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "partition_problem.hpp"

//****************************************************************************80

int i4_abs ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_ABS returns the absolute value of an I4.
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
//    Input, int I, an integer.
//
//    Output, int I4_ABS, the absolute value of the integer.
//
{
  int value;

  if ( 0 <= i )
  {
    value = i;
  }
  else
  {
    value = - i;
  }
  return value;
}
//****************************************************************************80

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

int *i4vec_copy_new ( int n, int a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY_NEW copies an I4VEC to a "new" I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int I4VEC_COPY_NEW[N], the copy of A1.
//
{
  int *a2;
  int i;

  a2 = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

int i4vec_dot_product ( int n, int x[], int y[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, int X[N], Y[N], the arrays.
//
//    Output, int I4VEC_DOT_PRODUCT, the dot product of X and Y.
//
{
  int i;
  int value;

  value = 0;
  for ( i = 0; i < n; i++ )
  {
    value = value + x[i] * y[i];
  }

  return value;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
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
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

void partition_brute ( int n, int w[], int c[], int &discrepancy )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_BRUTE approaches the partition problem using brute force.
//
//  Discussion:
//
//    We are given a set of N integers W.
//
//    We seek to partition W into subsets W0 and W1, such that the subsets
//    have equal sums.
//
//    The "discrepancy" is the absolute value of the difference between the
//    two sums, and will be zero if we have solved the problem.
//
//    For a given set of integers, there may be zero, one, or many solutions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the set.
//
//    Input, int W[N], the integers.
//
//    Output, int C[N], indicates the proposed solution.
//    C(I) is 0 for items in set W0 and 1 for items in set W1.
//
//    Output, int &DISCREPANCY, the discrepancy.
//
{
  int *d;
  int d_discrepancy;
  int rank;
  int w_sum;

  w_sum = i4vec_sum ( n, w );
  discrepancy = w_sum;

  rank = -1;
  d = new int[n];

  while ( 1 )
  {
    subset_next ( n, d, rank );

    if ( rank == -1 )
    {
      break;
    }

    d_discrepancy = i4_abs ( w_sum - 2 * i4vec_dot_product ( n, d, w ) );

    if ( d_discrepancy < discrepancy )
    {
      discrepancy = d_discrepancy;
      i4vec_copy ( n, d, c );
    }

    if ( discrepancy == 0 )
    {
      break;
    }
  }
  delete [] d;

  return;
}
//****************************************************************************80

int partition_count ( int n, int w[] )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_COUNT counts the solutions to a partition problem.
//
//  Discussion:
//
//    We are given a set of N integers W.
//
//    We seek to partition W into subsets W0 and W1, such that the subsets
//    have equal sums.
//
//    The "discrepancy" is the absolute value of the difference between the
//    two sums, and will be zero if we have solved the problem.
//
//    For a given set of integers, there may be zero, one, or many solutions.
//
//    In the case where the weights are distinct, the count returned by this
//    function may be regarded as twice as big as it should be, since the
//    partition (W0,W1) is counted a second time as (W1,W0).  A more serious
//    overcount can occur if the set W contains duplicate elements - in the
//    extreme case, W might be entirely 1's, in which case there is really
//    only one (interesting) solution, but this function will count many.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the set.
//
//    Input, int W[N], the integers.
//
//    Output, int PARTITION_COUNT, the number of solutions.
//
{
  int *c;
  int count;
  int discrepancy;
  int rank;
  int w_sum;

  w_sum = i4vec_sum ( n, w );

  c = new int[n];
  rank = -1;
  count = 0;

  while ( 1 )
  {
    subset_next ( n, c, rank );

    if ( rank == -1 )
    {
      break;
    }

    discrepancy = i4_abs ( w_sum - 2 * i4vec_dot_product ( n, c, w ) );

    if ( discrepancy == 0 )
    {
      count = count + 1;
    }
  }

  delete [] c;

  return count;
}
//****************************************************************************80

void subset_next ( int n, int t[], int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_NEXT computes the subset lexicographic successor.
//
//  Discussion:
//
//    This is a lightly modified version of "subset_lex_successor()" from COMBO.
//
//  Example:
//
//    On initial call, N is 5 and the input value of RANK is -1.
//    Then here are the successive outputs from the program:
//
//   Rank   T1   T2   T3   T4   T5
//   ----   --   --   --   --   --
//      0    0    0    0    0    0
//      1    0    0    0    0    1
//      2    0    0    0    1    0
//      3    0    0    0    1    1
//     ..   ..   ..   ..   ..   ..
//     30    1    1    1    1    0
//     31    1    1    1    1    1
//     -1    0    0    0    0    0  <-- Reached end of cycle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
//
//  Parameters:
//
//    Input, int N, the number of elements in the master set.
//    N must be positive.
//
//    Input/output, int T[N], describes a subset.  T(I) is 0 if
//    the I-th element of the master set is not in the subset, and is
//    1 if the I-th element is part of the subset.
//    On input, T describes a subset.
//    On output, T describes the next subset in the ordering.
//
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is -1.
//
{
  int i;
//
//  Return the first element.
//
  if ( rank == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      t[i] = 0;
    }
    rank = 0;
    return;
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( t[i] == 0 )
    {
      t[i] = 1;
      rank = rank + 1;
      return;
    }
    else
    {
      t[i] = 0;
    }
  }
  rank = -1;

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
