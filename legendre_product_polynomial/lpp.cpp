# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lpp.hpp"

//****************************************************************************80

int comp_enum ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_ENUM returns the number of compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The 28 compositions of 6 into three parts are:
//
//      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
//      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
//      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
//      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
//      0 3 3,  0 2 4,  0 1 5,  0 0 6.
//
//    The formula for the number of compositions of N into K parts is
//
//      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
//
//    Describe the composition using N '1's and K-1 dividing lines '|'.
//    The number of distinct permutations of these symbols is the number
//    of compositions.  This is equal to the number of permutations of
//    N+K-1 things, with N identical of one kind and K-1 identical of another.
//
//    Thus, for the above example, we have:
//
//      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2013
//
//  Author:
//
//    John Burkardt
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
//    Output, int COMP_ENUM, the number of compositions of N
//    into K parts.
//
{
  int number;

  number = i4_choose ( n + k - 1, n );

  return number;
}
//****************************************************************************80

void comp_next_grlex ( int kc, int xc[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT_GRLEX returns the next composition in grlex order.
//
//  Discussion:
//
//    Example:
//
//    KC = 3
//
//    #   XC(1) XC(2) XC(3)  Degree
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
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int KC, the number of parts of the composition.
//    1 <= KC.
//
//    Input/output, int XC[KC], the current composition.
//    Each entry of XC must be nonnegative.
//    On return, XC has been replaced by the next composition in the
//    grlex order.
//
{
  int i;
  int im1;
  int j;
  int t;
//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_NEXT_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= XC(I).
//
  for ( i = 0; i < kc; i++ )
  {
    if ( xc[i] < 0 )
    {
      cerr << "\n";
      cerr << "COMP_NEXT_GRLEX - Fatal error!\n";
      cerr << "  XC[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  Find I, the index of the rightmost nonzero entry of X.
//
  i = 0;
  for ( j = kc; 1 <= j; j-- )
  {
    if ( 0 < xc[j-1] )
    {
      i = j;
      break;
    }
  }
//
//  set T = X(I)
//  set XC(I) to zero,
//  increase XC(I-1) by 1,
//  increment XC(KC) by T-1.
//
  if ( i == 0 )
  {
    xc[kc-1] = 1;
    return;
  }
  else if ( i == 1 )
  {
    t = xc[0] + 1;
    im1 = kc;
  }
  else if ( 1 < i )
  {
    t = xc[i-1];
    im1 = i - 1;
  }

  xc[i-1] = 0;
  xc[im1-1] = xc[im1-1] + 1;
  xc[kc-1] = xc[kc-1] + t - 1;

  return;
}
//****************************************************************************80

int *comp_random_grlex ( int kc, int rank1, int rank2, int &seed, int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANDOM_GRLEX: random composition with degree less than or equal to NC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int KC, the number of parts in the composition.
//
//    Input, int RANK1, RANK2, the minimum and maximum ranks.
//    1 <= RANK1 <= RANK2.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, int &RANK, the rank of the composition.
//
//    Output, int COMP_RANDOM_GRLEX[KC], the random composition.
//
{
  int *xc;

//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_RANDOM_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK1.
//
  if ( rank1 < 1 )
  {
    cerr << "\n";
    cerr << "COMP_RANDOM_GRLEX - Fatal error!\n";
    cerr << "  RANK1 < 1\n";
    exit ( 1 );
  }
//
//  Ensure that RANK1 <= RANK2.
//
  if ( rank2 < rank1 )
  {
    cerr << "\n";
    cerr << "COMP_RANDOM_GRLEX - Fatal error!\n";
    cerr << "  RANK2 < RANK1\n";
    exit ( 1 );
  }
//
//  Choose RANK between RANK1 and RANK2.
//
  rank = i4_uniform_ab ( rank1, rank2, seed );
//
//  Recover the composition of given RANK.
//
  xc = comp_unrank_grlex ( kc, rank );

  return xc;
}
//****************************************************************************80

int comp_rank_grlex ( int kc, int xc[] )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_RANK_GRLEX computes the graded lexicographic rank of a composition.
//
//  Discussion:
//
//    The graded lexicographic ordering is used, over all KC-compositions
//    for NC = 0, 1, 2, ...
//
//    For example, if KC = 3, the ranking begins:
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
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int KC, the number of parts in the composition.
//    1 <= KC.
//
//    Input, int XC[KC], the composition.
//    For each 1 <= I <= KC, we have 0 <= XC(I).
//
//    Output, int COMP_RANK_GRLEX, the rank of the composition.
//
{
  int i;
  int j;
  int ks;
  int n;
  int nc;
  int ns;
  int rank;
  int tim1;
  int *xs;
//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_RANK_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 0 <= XC(I).
//
  for ( i = 0; i < kc; i++ )
  {
    if ( xc[i] < 0 )
    {
      cerr << "\n";
      cerr << "COMP_RANK_GRLEX - Fatal error!\n";
      cerr << "  XC[I] < 0\n";
      exit ( 1 );
    }
  }
//
//  NC = sum ( XC )
//
  nc = i4vec_sum ( kc, xc );
//
//  Convert to KSUBSET format.
//
  ns = nc + kc - 1;
  ks = kc - 1;
  xs = new int[ks];
  xs[0] = xc[0] + 1;
  for ( i = 2; i < kc; i++ )
  {
    xs[i-1] = xs[i-2] + xc[i-1] + 1;
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

  for ( n = 0; n < nc; n++ )
  {
    rank = rank + i4_choose ( n + kc - 1, n );
  }

  delete [] xs;

  return rank;
}
//****************************************************************************80

int *comp_unrank_grlex ( int kc, int rank )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_UNRANK_GRLEX computes the composition of given grlex rank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int KC, the number of parts of the composition.
//    1 <= KC.
//
//    Input, int RANK, the rank of the composition.
//    1 <= RANK.
//
//    Output, int COMP_UNRANK_GRLEX[KC], the composition XC of the given rank.
//    For each I, 0 <= XC[I] <= NC, and 
//    sum ( 1 <= I <= KC ) XC[I] = NC.
//
{
  int i;
  int j;
  int ks;
  int nc;
  int nksub;
  int ns;
  int r;
  int rank1;
  int rank2;
  int *xc;
  int *xs;
//
//  Ensure that 1 <= KC.
//
  if ( kc < 1 )
  {
    cerr << "\n";
    cerr << "COMP_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  KC < 1\n";
    exit ( 1 );
  }
//
//  Ensure that 1 <= RANK.
//
  if ( rank < 1 )
  {
    cerr << "\n";
    cerr << "COMP_UNRANK_GRLEX - Fatal error!\n";
    cerr << "  RANK < 1\n";
    exit ( 1 );
  }
//
//  Determine the appropriate value of NC.
//  Do this by adding up the number of compositions of sum 0, 1, 2, 
//  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
//  gives you the rank of the composition within the set of compositions
//  of sum NC.  And that's the number you need in order to do the
//  unranking.
//
  rank1 = 1;
  nc = -1;
  for ( ; ; )
  {
    nc = nc + 1;
    r = i4_choose ( nc + kc - 1, nc );
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
  ks = kc - 1;
  ns = nc + kc - 1;
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
  xc = new int[kc];
  xc[0] = xs[0] - 1;
  for ( i = 2; i < kc; i++ )
  {
    xc[i-1] = xs[i-1] - xs[i-2] - 1;
  }
  xc[kc-1] = ns - xs[ks-1];

  delete [] xs;

  return xc;
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
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
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
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
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

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
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

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
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
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
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
    cout << "  " << setw(8) << i
         << ": " << setw(8) << a[i]  << "\n";
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
//      i4vec_permute ( n, indx, a )
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
//    29 May 2003
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

int *i4vec_uniform_ab_new ( int n, int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_AB_NEW returns a scaled pseudorandom I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The pseudorandom numbers should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 May 2012
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
//    Input, int N, the dimension of the vector.
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int IVEC_UNIFORM_AB_NEW[N], a vector of random values 
//    between A and B.
//
{
  int c;
  int i;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;
  int *x;
  
  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  x = new int[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
      +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
    value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
    if ( value < a )
    {
      value = a;
    }
    if ( b < value )
    {
      value = b;
    }

    x[i] = value;
  }

  return x;
}
//****************************************************************************80

void lp_coefficients ( int n, int &o, double c[], int f[] )

//****************************************************************************80
//
//  Purpose:
//
//    LP_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
//
//  First terms:
//
//     1
//     0     1
//    -1/2   0      3/2
//     0    -3/2    0     5/2
//     3/8   0    -30/8   0     35/8
//     0    15/8    0   -70/8    0     63/8
//    -5/16  0    105/16  0   -315/16   0    231/16
//     0   -35/16   0   315/16   0   -693/16   0    429/16
//
//     1.00000
//     0.00000  1.00000
//    -0.50000  0.00000  1.50000
//     0.00000 -1.50000  0.00000  2.5000
//     0.37500  0.00000 -3.75000  0.00000  4.37500
//     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
//    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
//     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2014
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
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the highest order polynomial to evaluate.
//    Note that polynomials 0 through N will be evaluated.
//
//    Output, int &O, the number of coefficients.
//
//    Output, double C[(N+2)/2], the coefficients of the Legendre
//    polynomial of degree N.
//
//    Output, int F[(N+2)/2], the exponents.
//
{
  double *ctable;
  int i;
  int j;
  int k;

  ctable = new double[(n+1)*(n+1)];

  for ( i = 0; i <= n; i++ )
  {
    for ( j = 0; j <= n; j++ )
    {
      ctable[i+j*(n+1)] = 0.0;
    }
  }

  ctable[0+0*(n+1)] = 1.0;

  if ( 0 < n )
  {
    ctable[1+1*(n+1)] = 1.0;

    for ( i = 2; i <= n; i++ )
    {
      for ( j = 0; j <= i-2; j++ )
      {
        ctable[i+j*(n+1)] =
            ( double ) (   - i + 1 ) * ctable[i-2+j*(n+1)] / ( double ) i;
      }
      for ( j = 1; j <= i; j++ )
      {
        ctable[i+j*(n+1)] = ctable[i+j*(n+1)]
          + ( double ) ( i + i - 1 ) * ctable[i-1+(j-1)*(n+1)] / ( double ) i;
      }
    }
  }
//
//  Extract the nonzero data from the alternating columns of the last row.
//
  o = ( n + 2 ) / 2;

  k = o;
  for ( j = n; 0 <= j; j = j - 2 )
  {
    k = k - 1;
    c[k] = ctable[n+j*(n+1)];
    f[k] = j;
  }

  delete [] ctable;

  return;
}
//****************************************************************************80

double *lp_value ( int n, int o, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LP_VALUE evaluates the Legendre polynomials P(n,x).
//
//  Discussion:
//
//    P(n,1) = 1.
//    P(n,-1) = (-1)^N.
//    | P(n,x) | <= 1 in [-1,1].
//
//    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
//    quadrature of the integral of a function F(X) with weight function 1
//    over the interval [-1,1].
//
//    The Legendre polynomials are orthogonal under the inner product defined
//    as integration from -1 to 1:
//
//      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
//        = 0 if I =/= J
//        = 2 / ( 2*I+1 ) if I = J.
//
//    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
//
//    A function F(X) defined on [-1,1] may be approximated by the series
//      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
//    where
//      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
//
//    The formula is:
//
//      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
//
//  Differential equation:
//
//    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
//
//  First terms:
//
//    P( 0,x) =      1
//    P( 1,x) =      1 X
//    P( 2,x) = (    3 X^2 -       1)/2
//    P( 3,x) = (    5 X^3 -     3 X)/2
//    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
//    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
//    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
//    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
//    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
//    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
//    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
//
//  Recursion:
//
//    P(0,x) = 1
//    P(1,x) = x
//    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
//
//    P'(0,x) = 0
//    P'(1,x) = 1
//    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2014
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
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, int O, the degree of the polynomial.
//
//    Input, double X[N], the evaluation points.
//
//    Output, double LP_VALUE[N], the value of the Legendre polynomial 
//    of degree N at the points X.
//
{
  int i;
  int j;
  double *v;
  double *vtable;

  vtable = new double[n*(o+1)];

  for ( i = 0; i < n; i++ )
  {
    vtable[i+0*n] = 1.0;
  }

  if ( 1 <= o )
  {
    for ( i = 0; i < n; i++ )
    {
      vtable[i+1*n] = x[i];
    }

    for ( j = 2; j <= o; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        vtable[i+j*n] = 
          ( ( double ) ( 2 * j - 1 ) * x[i] * vtable[i+(j-1)*n]
          - ( double ) (     j - 1 ) *        vtable[i+(j-2)*n] )
          / ( double ) (     j     );
      }
    }
  }

  v = new double[n];

  for ( i = 0; i < n; i++ )
  {
    v[i] = vtable[i+o*n];
  }

  delete [] vtable;

  return v;
}
//****************************************************************************80

void lp_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LP_VALUES returns values of the Legendre polynomials P(n,x).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2014
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
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.4062500000000000E+00,
     -0.3359375000000000E+00,
      0.1577148437500000E+00,
      0.3397216796875000E+00,
      0.2427673339843750E-01,
     -0.2799186706542969E+00,
     -0.1524540185928345E+00,
      0.1768244206905365E+00,
      0.2212002165615559E+00,
      0.0000000000000000E+00,
     -0.1475000000000000E+00,
     -0.2800000000000000E+00,
     -0.3825000000000000E+00,
     -0.4400000000000000E+00,
     -0.4375000000000000E+00,
     -0.3600000000000000E+00,
     -0.1925000000000000E+00,
      0.8000000000000000E-01,
      0.4725000000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3,
     3 };

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
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     1.00E+00 };

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

void lpp_to_polynomial ( int m, int l[], int o_max, int &o, double c[], int e[] )

//****************************************************************************80
//
//  Purpose:
//
//    LPP_TO_POLYNOMIAL writes a Legendre Product Polynomial as a polynomial.
//
//  Discussion:
//
//    For example, if 
//      M = 3,
//      L = ( 1, 0, 2 ),
//    then
//      L(1,0,2)(X,Y,Z) 
//      = L(1)(X) * L(0)(Y) * L(2)(Z)
//      = X * 1 * ( 3Z^2-1)/2
//      = - 1/2 X + (3/2) X Z^2
//    so
//      O = 2 (2 nonzero terms)
//      C = -0.5
//           1.5
//      E = 4    <-- index in 3-space of exponent (1,0,0)
//          15   <-- index in 3-space of exponent (1,0,2)
//
//    The output value of O is no greater than
//      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int L[M], the index of each Legendre product 
//    polynomial factor.  0 <= L(*).
//
//    Input, int O_MAX, an upper limit on the size of the 
//    output arrays.
//      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
//
//    Output, int &O, the "order" of the polynomial product.
//
//    Output, double C[O], the coefficients of the polynomial product.
//
//    Output, int E[O], the indices of the exponents of the 
//    polynomial product.
//
{
  double *c1;
  double *c2;
  int *e1;
  int *e2;
  int *f2;
  int i;
  int i1;
  int i2;
  int j1;
  int j2;
  int o1;
  int o2;
  int *p;
  int *pp;

  c1 = new double[o_max];
  c2 = new double[o_max];
  e1 = new int[o_max];
  e2 = new int[o_max];
  f2 = new int[o_max];
  pp = new int[m];

  o1 = 1;
  c1[0] = 1.0;
  e1[0] = 1;
//
//  Implicate one factor at a time.
//
  for ( i = 0; i < m; i++ )
  {
    lp_coefficients ( l[i], o2, c2, f2 );
 
    o = 0;

    for ( j2 = 0; j2 < o2; j2++ )
    {
      for ( j1 = 0; j1 < o1; j1++ )
      {
        c[o] = c1[j1] * c2[j2];
        if ( 0 < i )
        {
          p = mono_unrank_grlex ( i, e1[j1] );
        }
        for ( i2 = 0; i2 < i; i2++ )
        {
          pp[i2] = p[i2];
        }
        pp[i] = f2[j2];
        e[o] = mono_rank_grlex ( i + 1, pp );
        o = o + 1;
        if ( 0 < i )
        {
          delete [] p;
        }
      }
    }

    polynomial_sort ( o, c, e );
    polynomial_compress ( o, c, e, o, c, e );

    o1 = o;
    for ( i1 = 0; i1 < o; i1++ )
    {
      c1[i1] = c[i1];
      e1[i1] = e[i1];
    }
  }

  delete [] c1;
  delete [] c2;
  delete [] e1;
  delete [] e2;
  delete [] f2;
  delete [] pp;

  return;
}
//****************************************************************************80

double *lpp_value ( int m, int n, int o[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LPP_VALUE evaluates a Legendre Product Polynomial at several points X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2014
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
//    Input, int O[M], the degree of the polynomial factors.
//    0 <= O(*).
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double LPP_VALUE[N], the value of the Legendre Product 
//    Polynomial of degree O at the points X.
//
{
  int i;
  int j;
  double *v;
  double *vi;
  double *xi;

  v = new double[n];

  for ( j = 0; j < n; j++ )
  {
    v[j] = 1.0;
  }

  xi = new double[n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      xi[j] = x[i+j*m];
    }
    vi = lp_value ( n, o[i], xi );
    for ( j = 0; j < n; j++ )
    {
      v[j] = v[j] * vi[j];
    }
    delete [] vi;
  }

  delete [] xi;

  return v;
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

void mono_print ( int m, int f[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_PRINT prints a monomial.
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
//    Input, int F[M], the exponents.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << title;
  cout << "x^(";
  for ( i = 0; i < m; i++ )
  {
    cout << f[i];
    if ( i < m - 1 )
    {
      cout << ",";
    }
    else
    {
      cout << ").\n";
    }
  }

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

void mono_upto_next_grlex ( int m, int n, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_NEXT_GRLEX: grlex next monomial with total degree up to N.
//
//  Discussion:
//
//    We consider all monomials in a M dimensional space, with total
//    degree up to N.
//
//    For example:
//
//    M = 3
//    N = 3
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
//    08 December 2013
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
//    0 <= N.
//
//    Input/output, int X[M], the current monomial.
//    To start the sequence, set X = [ 0, 0, ..., 0, 0 ].
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
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  if ( i4vec_sum ( m, x ) < 0 )
  {
    cerr << "\n";
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to less than 0.\n";
    exit ( 1 );
  }

  if ( n < i4vec_sum ( m, x ) )
  {
    cerr << "\n";
    cerr << "MONO_UPTO_NEXT_GRLEX - Fatal error!\n";
    cerr << "  Input X sums to more than N.\n";
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

int *mono_upto_random ( int m, int n, int &seed, int &rank )

//****************************************************************************80
//
//  Purpose:
//
//    MONO_UPTO_RANDOM: random monomial with total degree less than or equal to N.
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
//    Input, int M, the spatial dimension.
//
//    Input, int N, the degree.
//    0 <= N.
//
//    Input/output, int &SEED, the random number seed.
//
//    Output, int &RANK, the rank of the monomial.
//
//    Output, int MONO_UPTO_RANDOM[M], the random monomial.
//
{
  int rank_max;
  int rank_min;
  int *x;

  rank_min = 1;
  rank_max = mono_upto_enum ( m, n );
  rank = i4_uniform_ab ( rank_min, rank_max, seed );
  x = mono_unrank_grlex ( m, rank );

  return x;
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

int *perm_uniform_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of
//    (0, 1, ..., N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
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
//    The function polynomial_sort ( ) should be called first.
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

  get = 0;
  put = 0;

  while ( get < o1 )
  {
    get = get + 1;

    if ( fabs ( c1[get-1] ) <= r8_epsilon_sqrt )
    {
      continue;
    }

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
  }

  return p;
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
//    26 June 2013
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
    if ( n < j2hi )
    {
      j2hi = n;
    }
    if ( jhi < j2hi )
    {
      j2hi = jhi;
    }
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
    if ( 1 < ilo )
    {
      i2lo = ilo;
    }
    else
    {
      i2lo = 1;
    }
    if ( ihi < m )
    {
      i2hi = ihi;
    }
    else
    {
      i2hi = m;
    }

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

double *r8mat_uniform_ab_new ( int m, int n, double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB_NEW returns a new scaled pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    09 April 2012
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
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

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 May 2003
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
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << "  " << setw(8)  << i 
         << "  " << setw(12) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
//
//  Discussion:
//
//    Each dimension ranges from A to B.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    09 April 2012
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
//    Input, int N, the number of entries in the vector.
//
//    Input, double A, B, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
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
