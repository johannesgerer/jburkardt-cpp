# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>
# include <ctime>

using namespace std;

# include "polynomial.hpp"

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
//    06 January 2013
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

  mn = k;
  if ( n - k < mn )
  {
    mn = n - k;
  }

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
    mx = k;
    if ( mx < n - k )
    {
      mx = n - k;
    }
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

int i4_fall ( int x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FALL computes the falling factorial function [X]_N.
//
//  Discussion:
//
//    Note that the number of "injections" or 1-to-1 mappings from
//    a set of N elements to a set of M elements is [M]_N.
//
//    The number of permutations of N objects out of M is [M]_N.
//
//    Moreover, the Stirling numbers of the first kind can be used
//    to convert a falling factorial into a polynomial, as follows:
//
//      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
//
//    The formula is:
//
//      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X, the argument of the falling factorial function.
//
//    Input, int N, the order of the falling factorial function.
//    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
//    negative, a "rising" factorial will be computed.
//
//    Output, int I4_FALL, the value of the falling factorial function.
//
{
  int i;
  int value;

  value = 1;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * x;
      x = x - 1;
    }
  }
  else if ( n < 0 )
  {
    for ( i = -1; n <= i; i-- )
    {
      value = value * x;
      x = x + 1;
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

void i4vec_concatenate ( int n1, int a[], int n2, int b[], int c[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_CONCATENATE concatenates two I4VEC's.
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
//    22 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, the number of entries in the first vector.
//
//    Input, int A[N1], the first vector.
//
//    Input, int N2, the number of entries in the second vector.
//
//    Input, int B[N2], the second vector.
//
//    Output, int C[N1+N2], the concatenated vector.
//
{
  int i;

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return;
}
//****************************************************************************80

void i4vec_permute ( int n, int p[], int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PERMUTE permutes an I4VEC in place.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    This routine permutes an array of integer "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = (   1,   2,   3,   4,   5 )
//
//    Output:
//
//      A    = (   2,   4,   5,   1,   3 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input/output, int A[N], the array to be permuted.
//
{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  perm_check0 ( n, p );
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is 0.
//  So temporarily add 1 to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "I4VEC_PERMUTE - Fatal error!\n";
          cerr << "  Entry IPUT = " << iput << " of the permutation has\n";
          cerr << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
  }

  return;
}
//****************************************************************************80

int *i4vec_sort_heap_index_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(indx(*))
//
//    or explicitly, by the call
//
//      i4vec_permute ( n, indx, 0, a )
//
//    after which a(*) is sorted.
//
//    Note that the index vector is 0-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], an array to be index-sorted.
//
//    Output, int I4VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  int aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0];
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {

    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
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

void mono_next_grlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_NEXT_GRLEX returns the next monomial in grlex order.
//
//  Discussion:
//
//    Example:
//
//    M = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     0        0
//      |
//    2 |  0     0     1        1
//    3 |  0     1     0        1
//    4 |  1     0     0        1
//      |
//    5 |  0     0     2        2
//    6 |  0     1     1        2
//    7 |  0     2     0        2
//    8 |  1     0     1        2
//    9 |  1     1     0        2
//   10 |  2     0     0        2
//      |
//   11 |  0     0     3        3
//   12 |  0     1     2        3
//   13 |  0     2     1        3
//   14 |  0     3     0        3
//   15 |  1     0     2        3
//   16 |  1     1     1        3
//   17 |  1     2     0        3
//   18 |  2     0     1        3
//   19 |  2     1     0        3
//   20 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int X[M], the current monomial.
//    The first element is X = [ 0, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_NEXT_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  Find I, the index of the rightmost nonzero entry of X.
//
  i = 0;
  for ( j = m; 1 <= j; j-- )
  {
    if ( 0 < x[j-1] )
    {
      i = j;
      break;
    }
  }
//
//  set T = X(I)
//  set X(I) to zero,
//  increase X(I-1) by 1,
//  increment X(D) by T-1.
//
  if ( i == 0 )
  {
    x[m-1] = 1;
    return;
  }
  else if ( i == 1 )
  {
    t = x[0] + 1;
    im1 = m;
  }
  else if ( 1 < i )
  {
    t = x[i-1];
    im1 = i - 1;
  }

  x[i-1] = 0;
  x[im1-1] = x[im1-1] + 1;
  x[m-1] = x[m-1] + t - 1;

  return;
}
//****************************************************************************80

int mono_rank_grlex ( int m, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_RANK_GRLEX computes the graded lexicographic rank of a monomial.
//
//  Discussion:
//
//    The graded lexicographic ordering is used, over all monomials of
//    dimension M, for degree NM = 0, 1, 2, ...
//
//    For example, if M = 3, the ranking begins:
//
//    Rank  Sum    1  2  3
//    ----  ---   -- -- --
//       1    0    0  0  0
//
//       2    1    0  0  1
//       3    1    0  1  0
//       4    1    1  0  1
//
//       5    2    0  0  2
//       6    2    0  1  1
//       7    2    0  2  0
//       8    2    1  0  1
//       9    2    1  1  0
//      10    2    2  0  0
//
//      11    3    0  0  3
//      12    3    0  1  2
//      13    3    0  2  1
//      14    3    0  3  0
//      15    3    1  0  2
//      16    3    1  1  1
//      17    3    1  2  0
//      18    3    2  0  1
//      19    3    2  1  0
//      20    3    3  0  0
//
//      21    4    0  0  4
//      ..   ..   .. .. ..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= D.
//
//    Input, int XC[M], the monomial.
//    For each 1 <= I <= M, we have 0 <= XC(I).
//
//    Output, int MONO_RANK_GRLEX, the rank.
//
{
  int i;
  int j;
  int ks;
  int n;
  int nm;
  int ns;
  int rank;
  int tim1;
  int *xs;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_RANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= X(I).
//
  for ( i = 0; i < m; i++ )
  {
    if ( x[i] < 0 )
    {
      cerr << "\n";
      cerr << "MONO_RANK_GRLEX - Fatal error!\n";
      cerr << "  X[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  NM = sum ( X )
//
  nm = i4vec_sum ( m, x );
//
//  Convert to KSUBSET format.
//
  ns = nm + m - 1;
  ks = m - 1;
  xs = new int[ks];
  xs[0] = x[0] + 1;
  for ( i = 2; i < m; i++ )
  {
    xs[i-1] = xs[i-2] + x[i-1] + 1;
  }
//
//  Compute the rank.
//
  rank = 1;

  for ( i = 1; i <= ks; i++ )
  {
    if ( i == 1 )
    {
      tim1 = 0;
    }
    else
    {
      tim1 = xs[i-2];
    }

    if ( tim1 + 1 <= xs[i-1] - 1 )
    {
      for ( j = tim1 + 1; j <= xs[i-1] - 1; j++ )
      {
        rank = rank + i4_choose ( ns - j, ks - i );
      }
    }
  }

  for ( n = 0; n < nm; n++ )
  {
    rank = rank + i4_choose ( n + m - 1, n );
  }

  delete [] xs;

  return rank;
}
//****************************************************************************80

void mono_total_next_grlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_TOTAL_NEXT_GRLEX: grlex next monomial with total degree equal to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total degree N.
//
//    For example:
//
//    M = 3
//    N = 3
//
//    #  X(1)  X(2)  X(3)  Degree
//      +------------------------
//    1 |  0     0     3        3
//    2 |  0     1     2        3
//    3 |  0     2     1        3
//    4 |  0     3     0        3
//    5 |  1     0     2        3
//    6 |  1     1     1        3
//    7 |  1     2     0        3
//    8 |  2     0     1        3
//    9 |  2     1     0        3
//   10 |  3     0     0        3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, N ].
//    The last value in the sequence is X = [ N, 0, ..., 0, 0 ].
//
{
  int i;
  int im1;
  int j;
  int t;

  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) != n )
  {
    cerr << "\n";
    cerr << "MONO_TOTAL_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X does not sum to N.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    return;
  }

  if ( x[0] == n )
  {
    x[0] = 0;
    x[m-1] = n;
  }
  else
  {
    mono_next_grlex ( m, x );
  }

  return;
}
//****************************************************************************80

int *mono_unrank_grlex ( int m, int rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UNRANK_GRLEX computes the composition of given grlex rank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= M.
//
//    Input, int RANK, the rank.
//    1 <= RANK.
//
//    Output, int MONO_UNRANK_GRLEX[M], the monomial X of the given rank.
//    For each I, 0 <= XC[I] <= NM, and 
//    sum ( 1 <= I <= M ) XC[I] = NM.
//
{
  int i;
  int j;
  int ks;
  int nksub;
  int nm;
  int ns;
  int r;
  int rank1;
  int rank2;
  int *x;
  int *xs;
//
//  Ensure that 1 <= M.
//
  if ( m < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  M < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK.
//
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "MONO_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    exit ( 1 );
  }
//
//  Special case M == 1.
//
  if ( m == 1 )
  {
    x = new int[m];
    x[0] = rank - 1;
    return x;
  }
//
//  Determine the appropriate value of NM.
//  Do this by adding up the number of compositions of sum 0, 1, 2, 
//  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
//  gives you the rank of the composition within the set of compositions
//  of sum NM.  And that's the number you need in order to do the
//  unranking.
//
  rank1 = 1;
  nm = -1;
  for ( ; ; )
  {
    nm = nm + 1;
    r = i4_choose ( nm + m - 1, nm );
    if ( rank < rank1 + r )
    {
      break;
    }
    rank1 = rank1 + r;
  }

  rank2 = rank - rank1;
//
//  Convert to KSUBSET format.
//  Apology: an unranking algorithm was available for KSUBSETS,
//  but not immediately for compositions.  One day we will come back
//  and simplify all this.
//
  ks = m - 1;
  ns = nm + m - 1;
  xs = new int[ks];

  nksub = i4_choose ( ns, ks );

  j = 1;

  for ( i = 1; i <= ks; i++ )
  {
    r = i4_choose ( ns - j, ks - i );

    while ( r <= rank2 && 0 < r )
    {
      rank2 = rank2 - r;
      j = j + 1;
      r = i4_choose ( ns - j, ks - i );
    }
    xs[i-1] = j;
    j = j + 1;
  }
//
//  Convert from KSUBSET format to COMP format.
//
  x = new int[m];
  x[0] = xs[0] - 1;
  for ( i = 2; i < m; i++ )
  {
    x[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  x[m-1] = ns - xs[ks-1];

  delete [] xs;

  return x;
}
//****************************************************************************80

int mono_upto_enum ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_ENUM enumerates monomials in M dimensions of degree up to N.
//
//  Discussion:
//
//    For M = 2, we have the following values:
//
//    N  VALUE
//
//    0    1
//    1    3
//    2    6
//    3   10
//    4   15
//    5   21
//
//    In particular, VALUE(2,3) = 10 because we have the 10 monomials:
//
//      1, x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the maximum degree.
//
//    Output, int MONO_UPTO_ENUM, the number of monomials in
//    M variables, of total degree N or less.
//
{
  int value;

  value = i4_choose ( n + m, n );

  return value;
}
//****************************************************************************80

double *mono_value ( int m, int n, int f[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_VALUE evaluates a monomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
//    Input, int F[M], the exponents of the monomial.
//
//    Input, double X[M*N], the coordinates of the evaluation points.
//
//    Output, double MONO_VALUE[N], the value of the monomial at X.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
    for ( i = 0; i < m; i++ )
    {
      v[j] = v[j] * pow ( x[i+j*m], f[i] );
    }
  }

  return v;
}
//****************************************************************************80

void perm_check0 ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK0 checks a 0-based permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 0 to
//    to N-1 occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
{
  bool error;
  int location;
  int value;

  for ( value = 0; value < n; value++ )
  {
    error = true;

    for ( location = 0; location < n; location++ )
    {
      if ( p[location] == value )
      {
        error = false;
        break;
      }
    }

    if ( error )
    {
      cerr << "\n";
      cerr << "PERM_CHECK0 - Fatal error!\n";
      cerr << "  Permutation is missing value " << value << "\n";
      exit ( 1 );
    }

  }

  return;
}
//****************************************************************************80

void polynomial_add ( int o1, double c1[], int e1[], int o2, double c2[], 
  int e2[], int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_ADD adds two polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int O1, the "order" of polynomial 1.
//
//    Input, double C1[O1], the coefficients of polynomial 1.
//
//    Input, int E1[O1], the indices of the exponents of 
//    polynomial 1.
//
//    Input, int O2, the "order" of polynomial 2.
//
//    Input, double C2[O2], the coefficients of polynomial 2.
//
//    Input, int E2[O2], the indices of the exponents of 
//    polynomial 2.
//
//    Output, int &O, the "order" of the polynomial sum.
//
//    Output, double C[O], the coefficients of the polynomial sum.
//
//    Output, int E[O], the indices of the exponents of 
//    the polynomial sum.
//
{
  o = o1 + o2;
  r8vec_concatenate ( o1, c1, o2, c2, c );
  i4vec_concatenate ( o1, e1, o2, e2, e );

  polynomial_sort ( o, c, e );
  polynomial_compress ( o, c, e, o, c, e );

  return;
}
//****************************************************************************80

void polynomial_axpy ( double s, int o1, double c1[], int e1[], int o2, 
  double c2[], int e2[], int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_AXPY adds a multiple of one polynomial to another.
//
//  Discussion:
//
//    P(X) = S * P1(X) + P2(X)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the multiplier of polynomial 1.
//
//    Input, int O1, the "order" of polynomial 1.
//
//    Input, double C1[O1], the coefficients of polynomial 1.
//
//    Input, int E1[O1], the indices of the exponents of 
//    polynomial 1.
//
//    Input, int O2, the "order" of polynomial 2.
//
//    Input, double C2[O2], the coefficients of polynomial 2.
//
//    Input, int E2[O2], the indices of the exponents of 
//    polynomial 2.
//
//    Output, int &O, the "order" of the polynomial sum.
//
//    Output, double C[O], the coefficients of the polynomial sum.
//
//    Output, int E[O], the indices of the exponents of 
//    the polynomial sum.
//
{
  double *c3;
  int *e3;
  int i;
  int o3;
  double *sc1;

  o3 = o1 + o2;

  c3 = new double[o3];
  e3 = new int[o3];
  sc1 = new double[o1];

  for ( i = 0; i < o1; i++ )
  {
    sc1[i] = s * c1[i];
  }
  r8vec_concatenate ( o1, sc1, o2, c2, c3 );
  i4vec_concatenate ( o1, e1, o2, e2, e3 );

  polynomial_sort ( o3, c3, e3 );
  polynomial_compress ( o3, c3, e3, o, c, e );

  delete [] c3;
  delete [] e3;
  delete [] sc1;

  return;
}
//****************************************************************************80

void polynomial_compress ( int o1, double c1[], int e1[], int &o2, double c2[], 
  int e2[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_COMPRESS compresses a polynomial.
//
//  Discussion:
//
//    The function polynomial_sort ( ) should be called first, or else
//    the E1 vector should be in ascending sorted order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int O1, the "order" of the polynomial.
//
//    Input, double C1[O1], the coefficients of the polynomial.
//
//    Input, int E1[O1], the indices of the exponents of 
//    the polynomial.
//
//    Output, int &O2, the "order" of the polynomial.
//
//    Output, double C2[O2], the coefficients of the polynomial.
//
//    Output, int E2[O2], the indices of the exponents of 
//    the polynomial.
//
{
  int get;
  int put;
  const double r8_epsilon_sqrt = 0.1490116119384766E-07;
//
//  Add coefficients associated with the same exponent.
//
  get = 0;
  put = 0;

  while ( get < o1 )
  {
    get = get + 1;

    if ( 0 == put )
    {
      put = put + 1;
      c2[put-1] = c1[get-1];
      e2[put-1] = e1[get-1];
    }
    else
    {
      if ( e2[put-1] == e1[get-1] )
      {
        c2[put-1] = c2[put-1] + c1[get-1];
      }
      else
      {
        put = put + 1;
        c2[put-1] = c1[get-1];
        e2[put-1] = e1[get-1];
       }
    }
  }
 
  o2 = put;
//
//  Clear out zeros and tiny coefficients.
//
  get = 0;
  put = 0;

  while ( get < o2 )
  {
    if ( r8_epsilon_sqrt < fabs ( c2[get] ) )
    {
      c2[put] = c2[get];
      e2[put] = e2[get];
      put = put + 1;
    }
    get = get + 1;
  }

  o2 = put;
  return;
}
//****************************************************************************80

void polynomial_dif ( int m, int o1, double c1[], int e1[], int dif[], 
  int &o2, double c2[], int e2[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_DIF differentiates a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int O1, the "order" of polynomial 1.
//
//    Input, double C1[O1], the coefficients of polynomial 1.
//
//    Input, int E1[O1], the indices of the exponents of 
//    polynomial 1.
//
//    Input, int DIF[M], indicates the number of 
//    differentiations in each component.
//
//    Output, int &O2, the "order" of the polynomial derivative.
//
//    Output, double C2[O2], the coefficients of the polynomial 
//    derivative.
//
//    Output, int E2[O2], the indices of the exponents of the
//    polynomial derivative.
//
{
  int *f1;
  int i;
  int j;

  o2 = o1;
  for ( j = 0; j < o1; j++ )
  {
    c2[j] = c1[j];
  }

  for ( j = 0; j < o1; j++ )
  {
    f1 = mono_unrank_grlex ( m, e1[j] );
    for ( i = 0; i < m; i++ )
    {
      c2[j] = c2[j] * i4_fall ( f1[i], dif[i] );
      f1[i] = i4_max ( f1[i] - dif[i], 0 );
    }
    e2[j] = mono_rank_grlex ( m, f1 );
    delete [] f1;
  }

  polynomial_sort ( o2, c2, e2 );

  polynomial_compress ( o2, c2, e2, o2, c2, e2 );

  return;
}
//****************************************************************************80

void polynomial_mul ( int m, int o1, double c1[], int e1[], int o2, double c2[],
  int e2[], int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_MUL multiplies two polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int O1, the "order" of polynomial 1.
//
//    Input, double C1[O1], the coefficients of polynomial 1.
//
//    Input, int E1[O1], the indices of the exponents of 
//    polynomial 1.
//
//    Input, int O2, the "order" of polynomial 2.
//
//    Input, double C2[O2], the coefficients of polynomial 2.
//
//    Input, int E2[O2], the indices of the exponents of 
//    polynomial 2.
//
//    Output, int &O, the "order" of the polynomial product.
//
//    Output, double C[O], the coefficients of the polynomial product.
//
//    Output, int E[O], the indices of the exponents of the 
//    polynomial product.
//
{
  int *f;
  int *f1;
  int *f2;
  int i;
  int j;
  int k;

  f = new int[m];

  o = 0;
  for ( j = 0; j < o2; j++ )
  {
    for ( i = 0; i < o1; i++ )
    {
      c[o] = c1[i] * c2[j];
      f1 = mono_unrank_grlex ( m, e1[i] );
      f2 = mono_unrank_grlex ( m, e2[j] );
      for ( k = 0; k < m; k++ )
      {
        f[k] = f1[k] + f2[k];
      }
      e[o] = mono_rank_grlex ( m, f );
      delete [] f1;
      delete [] f2;
      o = o + 1;
    }
  }

  delete [] f;

  polynomial_sort ( o, c, e );
  polynomial_compress ( o, c, e, o, c, e );

  return;
}
//****************************************************************************80

void polynomial_print ( int m, int o, double c[], int e[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_PRINT prints a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int O, the "order" of the polynomial, that is,
//    simply the number of terms.
//
//    Input, double C[O], the coefficients.
//
//    Input, int E[O], the indices of the exponents.
//
//    Input, string TITLE, a title.
//
{
  int *f;
  int i;
  int j;

  cout << title << "\n";

  if ( o == 0 )
  {
    cout << "      0.\n";
  }
  else
  {
    for ( j = 0; j < o; j++ )
    {
      cout << "    ";
      if ( c[j] < 0.0 )
      {
        cout << "- ";
      }
      else
      {
        cout << "+ ";
      }
      cout << fabs ( c[j] ) << " * x^(";

      f = mono_unrank_grlex ( m, e[j] );
      for ( i = 0; i < m; i++ )
      {
        cout << f[i];
        if ( i < m - 1 )
        {
          cout << ",";
        }
        else
        {
          cout << ")";
        }
      }
      delete [] f;

      if ( j == o - 1 )
      {
        cout << ".";
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void polynomial_scale ( double s, int m, int o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SCALE scales a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double S, the scale factor.
//
//    Input, int M, the spatial dimension.
//
//    Input, int O, the "order" of the polynomial.
//
//    Input/output, double C[O], the coefficients of the polynomial.
//
//    Input, int E[O], the indices of the exponents of the
//    polynomial.
//
{
  int i;

  for ( i = 0; i < o; i++ )
  {
    c[i] = c[i] * s;
  }

  return;
}
//****************************************************************************80

void polynomial_sort ( int o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_SORT sorts the information in a polynomial.
//
//  Discussion
//
//    The coefficients C and exponents E are rearranged so that 
//    the elements of E are in ascending order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int O, the "order" of the polynomial.
//
//    Input/output, double C[O], the coefficients of the polynomial.
//
//    Input/output, int E[O], the indices of the exponents of 
//    the polynomial.
//
{
  int *indx;

  indx = i4vec_sort_heap_index_a ( o, e );

  i4vec_permute ( o, indx, e );
  r8vec_permute ( o, indx, c );

  delete [] indx;

  return;
}
//****************************************************************************80

double *polynomial_value ( int m, int o, double c[], int e[], int n, 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYNOMIAL_VALUE evaluates a polynomial.
//
//  Discussion:
//
//    The polynomial is evaluated term by term, and no attempt is made to
//    use an approach such as Horner's method to speed up the process.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int O, the "order" of the polynomial.
//
//    Input, double C[O], the coefficients of the polynomial.
//
//    Input, int E(O), the indices of the exponents 
//    of the polynomial.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the coordinates of the evaluation points.
//
//    Output, double POLYNOMIAL_VALUE[N], the value of the polynomial at X.
//
{
  int *f;
  int j;
  int k;
  double *p;
  double *v;

  p = new double[n];

  for ( k = 0; k < n; k++ )
  {
    p[k] = 0.0;
  }

  for ( j = 0; j < o; j++ )
  {
    f = mono_unrank_grlex ( m, e[j] );
    v = mono_value ( m, n, f, x );
    for ( k = 0; k < n; k++ )
    {
      p[k] = p[k] + c[j] * v[k];
    }
    delete [] f;
    delete [] v;
  }

  return p;
}
//****************************************************************************80

void r8vec_concatenate ( int n1, double a[], int n2, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CONCATENATE concatenates two R8VEC's.
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
//    22 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, the number of entries in the first vector.
//
//    Input, double A[N1], the first vector.
//
//    Input, int N2, the number of entries in the second vector.
//
//    Input, double B[N2], the second vector.
//
//    Output, double C[N1+N2], the concatenated vector.
//
{
  int i;

  for ( i = 0; i < n1; i++ )
  {
    c[i] = a[i];
  }
  for ( i = 0; i < n2; i++ )
  {
    c[n1+i] = b[i];
  }

  return;
}
//****************************************************************************80

void r8vec_permute ( int n, int p[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PERMUTE permutes an R8VEC in place.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   1,   3,   4,   0,   2 )
//      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
//
//    Output:
//
//      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.
//
//    Input/output, double A[N], the array to be permuted.
//
{
  double a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  perm_check0 ( n, p );
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is 0.
//  So temporarily add 1 to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cerr << "\n";
          cerr << "R8VEC_PERMUTE - Fatal error!\n";
          cerr << "  A permutation index is out of range.\n";
          cerr << "  P(" << iput << ") = " << iget << "\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1;
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
