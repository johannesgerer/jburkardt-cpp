# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "combo.hpp"

//****************************************************************************80

void backtrack ( int l, int iarray[], int &indx, int &k, int &nstack, 
  int stack[], int maxstack )

//****************************************************************************80
// 
//  Purpose:
//
//    BACKTRACK supervises a backtrack search.
// 
//  Discussion:
// 
//    The routine builds a vector, one element at a time, which is
//    required to satisfy some condition.
// 
//    At any time, the partial vector may be discovered to be
//    unsatisfactory, but the routine records information about where the
//    last arbitrary choice was made, so that the search can be
//    carried out efficiently, rather than starting out all over again.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int L, the length of the completed candidate vector.
// 
//    Input/output, int IARRAY[L], the candidate vector.
// 
//    Input/output, int &INDX.
//    On input, set INDX = 0 to start a search.
//    On output:
//    1, a complete output vector has been determined.
//    2, candidates are needed.
//    3, no more possible vectors exist.
// 
//    Input/output, int &K, the current length of the candidate
//    vector.
// 
//    Input/output, int &NSTACK, the current length of the stack.
// 
//    Input/output, int STACK[MAXSTACK], a list of candidates
//    for positions 1 through K.
// 
//    Input, int MAXSTACK, the maximum length of the stack.
// 
{
// 
//  If this is the first call, request a candidate for position 1.
// 
  if ( indx == 0 )
  {
    k = 1;
    nstack = 0;
    indx = 2;
    return;
  }
// 
//  Examine the stack.
// 
  for ( ; ; )
  {
    nstack = nstack - 1;
// 
//  If there are candidates for position K, take the first available
//  one off the stack, and increment K.
// 
//  This may cause K to reach the desired value of L, in which case
//  we need to signal the user that a complete set of candidates
//  is being returned.
// 
    if ( stack[nstack] != 0 )
    {
      iarray[k-1] = stack[nstack-1];
      stack[nstack-1] = stack[nstack] - 1;

      if ( k != l )
      {
        k = k + 1;
        indx = 2;
      }
      else
      {
        indx = 1;
      }
      break;
    }
// 
//  If there are no candidates for position K, then decrement K.
//  If K is still positive, repeat the examination of the stack.
// 
    else
    {
      k = k - 1;

      if ( k <= 0 )
      {
        indx = 3;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

int bal_seq_check ( int n, int t[] )

//****************************************************************************80
//
//  Purpose:
//
//    BAL_SEQ_CHECK checks a balanced sequence.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Input, int T[2*N], a balanced sequence.
// 
//    Output, int BAL_SEQ_CHECK, error flag.
//    0, no error.
//    -1, N is not positive.
//    I+1, the I-th entry of T is illegal.
//    2*N+1, there are not the same number of 1's as 0's.
// 
{
  int i;
  int ierror;
  int one_count;
  int zero_count;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  one_count = 0;
  zero_count = 0;

  for ( i = 0; i < 2 * n; i++ )
  {
    if ( t[i] == 0 )
    {
      zero_count = zero_count + 1;
    }
    else if ( t[i] == 1 )
    {
      one_count = one_count + 1;
    }
    else
    {
      ierror = i + 1;
      return ierror;
    }

    if ( zero_count < one_count )
    {
      ierror = 1;
      return ierror;
    }

  }

  if ( one_count != zero_count )
  {
    ierror = 2 * n + 1;
  }

  return ierror;
}
//****************************************************************************80

int bal_seq_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    BAL_SEQ_ENUM enumerates the balanced sequences.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be nonnegative.
// 
//    Output, int BAL_SEQ_ENUM, the number of balanced sequences.
// 
{
  int value;

  value = i4_choose ( 2 * n, n ) / ( n + 1 );

  return value;
}
//****************************************************************************80

int bal_seq_rank ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    BAL_SEQ_RANK ranks a balanced sequence.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Input, int T[2*N], a balanced sequence.
// 
//    Output, int BAL_SEQ_RANK, the rank of the balanced sequence.
// 
{
  int ierror;
  int mxy;
  int rank;
  int x;
  int y;
// 
//  Check.
// 
  ierror = bal_seq_check ( n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "BAL_SEQ_RANK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  y = 0;
  rank = 0;

  for ( x = 1; x <= 2 * n - 1; x++ )
  {
    if ( t[x-1] == 0 )
    {
      y = y + 1;
    }
    else
    {
      mxy = mountain ( n, x, y + 1 );
      rank = rank + mxy;
      y = y - 1;
    }
  }

  return rank;
}
//****************************************************************************80

void bal_seq_successor ( int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    BAL_SEQ_SUCCESSOR computes the lexical balanced sequence successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Input/output, int T[2*N], on input, a balanced sequence,
//    and on output, its lexical successor.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
  int j;
  int open;
  int open_index;
  int slot;
  int slot_index;
  int slot_ones;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      t[i] = 0;
    }
    for ( i = n; i < 2 * n; i++ )
    {
      t[i] = 1;
    }
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = bal_seq_check ( n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "BAL_SEQ_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  After the I-th 0 there is a 'slot' with the capacity to
//  hold between 0 and I ones.
// 
//  The first element of the sequence has all the 1''s cowering
//  behind the N-th 0.
// 
//  We seek to move a 1 to the left, and to do it lexically,
//  we will move a 1 to the rightmost slot that is under capacity.
// 
//  Find the slot.
// 
  slot = 0;
  slot_index = 0;
  slot_ones = 0;

  open = 0;
  open_index = 0;

  for ( i = 1; i <= 2 * n; i++ )
  {
    if ( t[i-1] == 0 )
    {
      if ( 0 < slot )
      {
        if ( slot_ones < slot )
        {
          open = slot;
          open_index = slot_index;
        }
      }
      slot = slot + 1;
      slot_index = i;
    }
    else
    {
      slot_ones = slot_ones + 1;
    }
  }
// 
//  If OPEN is not 0, then preserve the string up to the OPEN-th 0,
//  preserve the 1''s that follow, but then write a 1, then
//  all the remaining 0's and all the remaining 1's.
// 
  if ( open != 0 )
  {
    j = open_index + 1;

    while ( t[j-1] == 1 )
    {
      j = j + 1;
    }

    t[j-1] = 1;

    for ( i = open + 1; i <= n; i++ )
    {
      j = j + 1;
      t[j-1] = 0;
    }

    for ( i = j + 1; i <= 2 * n; i++ )
    {
      t[i-1] = 1;
    }
  }
// 
//  If OPEN is 0, the last element was input.
//  Return the first one.
// 
  else
  {
    for ( i = 0; i < n; i++ )
    {
      t[i] = 0;
    }
    for ( i = n; i < 2 * n; i++ )
    {
      t[i] = 1;
    }
    rank = 0;
    return;
  }
  rank = rank + 1;

  return;
}
//****************************************************************************80

int *bal_seq_to_tableau ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    BAL_SEQ_TO_TABLEAU converts a balanced sequence to a 2 by N tableau.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Input, int T[2*N], a balanced sequence.
// 
//    Output, int BAL_SEQ_TO_TABLEAU[2*N], a 2 by N tableau.
// 
{
  int c[2];
  int i;
  int ierror;
  int r;
  int *tab;
// 
//  Check.
// 
  ierror = bal_seq_check ( n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "BAL_SEQ_TO_TABLEAU - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  tab = new int[2*n];

  c[0] = 0;
  c[1] = 0;

  for ( i = 1; i <= 2 * n; i++ )
  {
    r = t[i-1] + 1;
    c[r-1] = c[r-1] + 1;
    tab[r-1+(c[r-1]-1)*2] = i;
  }

  return tab;
}
//****************************************************************************80

int *bal_seq_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    BAL_SEQ_UNRANK unranks a balanced sequence.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int RANK, the rank of the balanced sequence.
// 
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Output, int BAL_SEQ_UNRANK[2*N], a balanced sequence.
// 
{
  int low;
  int m;
  int nseq;
  int *t;
  int x;
  int y;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "BAL_SEQ_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  nseq = bal_seq_enum ( n );

  if ( rank < 0 || nseq < rank )
  {
    cout << "\n";
    cout << "BAL_SEQ_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  t = new int[2*n];

  y = 0;
  low = 0;

  for ( x = 0; x < 2 * n; x++ )
  {
    m = mountain ( n, x + 1, y + 1 );

    if ( rank <= low + m - 1 )
    {
      y = y + 1;
      t[x] = 0;
    }
    else
    {
      low = low + m;
      y = y - 1;
      t[x] = 1;
    }
  }
  return t;
}
//****************************************************************************80

int *bell_numbers ( int m )

//****************************************************************************80
// 
//  Purpose:
//
//    BELL_NUMBERS computes the Bell numbers.
// 
//  Discussion:
// 
//    There are B(M) restricted growth functions of length M.
// 
//    There are B(M) partitions of a set of M objects.
// 
//    B(M) is the sum of the Stirling numbers of the second kind,
//    S(M,N), for N = 0 to M.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int M, indicates how many Bell numbers are to
//    compute.  M must be nonnegative.
// 
//    Output, int BELL_NUMBERS[M+1], the first M+1 Bell numbers.
// 
{
  int *b;
  int i;
  int j;

  b = new int[m+1];

  b[0] = 1;
  for ( j = 1; j <= m; j++ )
  {
    b[j] = 0;
    for ( i = 0; i < j; i++ )
    {
      b[j] = b[j] + i4_choose ( j - 1, i ) * b[i];
    }
  }
  return b;
}
//****************************************************************************80

void bell_values ( int *n_data, int *n, int *c )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_VALUES returns some values of the Bell numbers.
//
//  Discussion:
//
//    The Bell number B(N) is the number of restricted growth functions on N.
//
//    Note that the Stirling numbers of the second kind, S^m_n, count the
//    number of partitions of N objects into M classes, and so it is
//    true that
//
//      B(N) = S^1_N + S^2_N + ... + S^N_N.
//
//    The Bell numbers were named for Eric Temple Bell.
//
//    In Mathematica, the function can be evaluated by
//
//      Sum[StirlingS2[n,m],{m,1,n}]
//
//  Definition:
//
//    The Bell number B(N) is defined as the number of partitions (of
//    any size) of a set of N distinguishable objects.
//
//    A partition of a set is a division of the objects of the set into
//    subsets.
//
//  Examples:
//
//    There are 15 partitions of a set of 4 objects:
//
//      (1234),
//      (123) (4),
//      (124) (3),
//      (12) (34),
//      (12) (3) (4),
//      (134) (2),
//      (13) (24),
//      (13) (2) (4),
//      (14) (23),
//      (1) (234),
//      (1) (23) (4),
//      (14) (2) (3),
//      (1) (24) (3),
//      (1) (2) (34),
//      (1) (2) (3) (4).
//
//    and so B(4) = 15.
//
//  First values:
//
//     N         B(N)
//     0           1
//     1           1
//     2           2
//     3           5
//     4          15
//     5          52
//     6         203
//     7         877
//     8        4140
//     9       21147
//    10      115975
//
//  Recursion:
//
//    B(I) = sum ( 1 <= J <=I ) i4_choose ( I-1, J-1 ) * B(I-J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2003
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
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int *N, the order of the Bell number.
//
//    Output, int *C, the value of the Bell number.
//
{
# define N_MAX 11

  static int c_vec[N_MAX] = {
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *n = 0;
    *c = 0;
  }
  else
  {
    *n = n_vec[*n_data-1];
    *c = c_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double combin ( int n, int k )

//****************************************************************************80
// 
//  Purpose:
//
//    COMBIN computes the combinatorial coefficient C(N,K).
// 
//  Discussion:
// 
//    Real arithmetic is used, and C(N,K) is computed directly, via
//    Gamma functions, rather than recursively.
// 
//    C(N,K) is the number of distinct combinations of K objects
//    chosen from a set of N distinct objects.  A combination is
//    like a set, in that order does not matter.
// 
//    C(N,K) = N! / ( (N-K)! * K! )
// 
//  Example:
// 
//    The number of combinations of 2 things chosen from 5 is 10.
// 
//    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
// 
//    The actual combinations may be represented as:
// 
//      (1,2), (1,3), (1,4), (1,5), (2,3),
//      (2,4), (2,5), (3,4), (3,5), (4,5).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the value of N.
// 
//    Input, int K, the value of K.
// 
//    Output, double COMBIN, the value of C(N,K)
// 
{
  double arg;
  double fack;
  double facn;
  double facnmk;
  double value;

  if ( n < 0 )
  {
    value = 0.0;
  }
  else if ( k == 0 )
  {
    value = 1.0;
  }
  else if ( k == 1 )
  {
    value = ( double ) ( n );
  }
  else if ( 1 < k && k < n - 1 )
  {
    arg = ( double ) ( n + 1 );
    facn = r8_gamma_log ( arg );

    arg = ( double ) ( k + 1 );
    fack = r8_gamma_log ( arg );

    arg = ( double ) ( n - k + 1 );
    facnmk = r8_gamma_log ( arg );

    value = r8_nint ( exp ( facn - fack - facnmk ) );
  }
  else if ( k == n - 1 )
  {
    value = ( double ) ( n );
  }
  else if ( k == n )
  {
    value = 1.0;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void cycle_check ( int n, int ncycle, int t[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    CYCLE_CHECK checks a permutation in cycle form.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int N, the number of items permuted.
//    N must be positive.
// 
//    Input, int NCYCLE, the number of cycles.
//    1 <= NCYCLE <= N.
// 
//    Input, int T[N], INDEX[NCYCLE], describes the permutation
//    as a collection of NCYCLE cycles.  The first cycle is
//    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
//
{
  int i;
  int ifind;
  int iseek;
// 
//  N must be at least 1.
// 
  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CYCLE_CHECK - Fatal error!\n";
    cerr << "  N must be at least 1.\n";
    exit ( 1 );
  }
// 
//  1 <= NCYCLE <= N.
// 
  if ( ncycle < 1 || n < ncycle )
  {
    cerr << "\n";
    cerr << "CYCLE_CHECK - Fatal error!\n";
    cerr << "  1 <= NCYCLE <= N is required.\n";
    exit ( 1 );
  }
// 
//  1 <= INDEX(I) <= N.
//
  for ( i = 0; i < ncycle; i++ )
  {
    if ( index[i] < 1 || n < index[i] )
    {
      cerr << "\n";
      cerr << "CYCLE_CHECK - Fatal error!\n";
      cerr << "  1 <= INDEX[I] <= N is required.\n";
      cerr << "  But index[" << i << "] = " << index[i] << "\n";
      exit ( 1 );
    }
  }
// 
//  The INDEX(I)''s sum to N.
// 
  if ( i4vec_sum ( ncycle, index ) != n )
  {
    cerr << "\n";
    cerr << "CYCLE_CHECK - Fatal error!\n";
    cerr << "  INDEX entries must sum to N.\n";
    exit ( 1 );
  }
// 
//  1 <= T(I) <= N.
// 
  for ( i = 0; i < n; i++ )
  {
    if ( t[i] < 1 || n < t[i] )
    {
      cerr << "\n";
      cerr << "CYCLE_CHECK - Fatal error!\n";
      cerr << "  1 <= T[I] <= N is required.\n";
      exit ( 1 );
    }
  }
// 
//  Verify that every value from 1 to N occurs in T.
// 
  for ( iseek = 1; iseek <= n; iseek++ )
  {
    ifind = -1;

    for ( i = 0; i < n; i++ )
    {
      if ( t[i] == iseek )
      {
        ifind = i + 1;
        break;
      }
    }

    if ( ifind == -1 )
    {
      cerr << "\n";
      cerr << "CYCLE_CHECK - Fatal error!\n";
      cerr << "  The value " << iseek << " does not occur in T.\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************80

int *cycle_to_perm ( int n, int ncycle, int t[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    CYCLE_TO_PERM converts a permutation from cycle to array form.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of items permuted.
//    N must be positive.
// 
//    Input, int NCYCLE, the number of cycles.
//    1 <= NCYCLE <= N.
// 
//    Input, int T[N], INDEX[NCYCLE], describes the permutation
//    as a collection of NCYCLE cycles.  The first cycle is
//    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
// 
//    Output, int CYCLE_TO_PERM[N], describes the permutation using a
//    single array.  For each index I, I -> P(I).
// 
{
  int i;
  int j;
  int jhi;
  int jlo;
  int *p;
// 
//  Check.
// 
  cycle_check ( n, ncycle, t, index );

  p = new int[n];

  jhi = 0;

  for ( i = 1; i <= ncycle; i++ )
  {
    jlo = jhi + 1;
    jhi = jhi + index[i-1];

    for ( j = jlo; j <= jhi; j++ )
    {
      if ( j < jhi )
      {
        p[t[j-1]-1] = t[j];
      }
      else
      {
        p[t[j-1]-1] = t[jlo-1];
      }
    }
  }

  return p;
}
//****************************************************************************80

int dist_enum ( int k, int m )

//****************************************************************************80
// 
//  Purpose:
//
//    DIST_ENUM returns the number of distributions of indistinguishable objects.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int K, the number of distinguishable "slots".
// 
//    Input, int M, the number of indistinguishable objects.
// 
//    Output, int DIST_ENUM, the number of distributions of M
//    indistinguishable objects about K distinguishable slots.
// 
{
  int value;

  value = i4_choose ( m + k - 1, m );

  return value;
}
//****************************************************************************80

void dist_next ( int k, int m, int q[], bool &more )

//****************************************************************************80
// 
//  Purpose:
//
//    DIST_NEXT returns the next distribution of indistinguishable objects.
// 
//  Discussion:
// 
//    A distribution of M objects into K parts is an ordered sequence
//    of K nonnegative integers which sum to M.  This is similar to
//    a partition of a set into K subsets, except that here the order
//    matters.  That is, (1,1,2) and (1,2,1) are considered to be
//    different distributions.
// 
//    On the first call to this routine, the user should set MORE = FALSE,
//    to signal that this is a startup for the given computation.  The routine
//    will return the first distribution, and set MORE = TRUE.
// 
//    If the user calls again, with MORE = TRUE, the next distribution
//    is being requested.  If the routine returns with MORE = TRUE, then
//    that distribution was found and returned.  However, if the routine
//    returns with MORE = FALSE, then no more distributions were found;
//    the enumeration of distributions has terminated.
// 
//    A "distribution of M indistinguishable objects into K slots" is
//    sometimes called a "composition of the integer M into K parts".
// 
//  Example:
// 
//    K = 3, M = 5
// 
//    0           0           5
//    0           1           4
//    0           2           3
//    0           3           2
//    0           4           1
//    0           5           0
//    1           0           4
//    1           1           3
//    1           2           2
//    1           3           1
//    1           4           0
//    2           0           3
//    2           1           2
//    2           2           1
//    2           3           0
//    3           0           2
//    3           1           1
//    3           2           0
//    4           0           1
//    4           1           0
//    5           0           0
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Reference:
// 
//    Robert Fenichel,
//    Algorithm 329:
//    Distribution of Indistinguishable Objects into
//    Distinguishable Slots,
//    Communications of the ACM,
//    Volume 11, Number 6, June 1968, page 430.
// 
//  Parameters:
// 
//    Input, int K, the number of distinguishable "slots".
// 
//    Input, int M, the number of indistinguishable objects.
// 
//    Input/output, int Q[K], the number of objects in each
//    slot.
// 
//    Input/output, bool &MORE, used by the user to start the computation,
//    and by the routine to stop the computation.
// 
{
  int i;
  static int leftmost = 1;
// 
//  The startup call.
// 
  if ( !more )
  {
    more = true;
    for ( i = 0; i < k - 1; i++ )
    {
      q[i] = 0;
    }
    q[k-1] = m;

    leftmost = k + 1;
  }
// 
//  There are no more distributions.
//  Reset Q to the first distribution in the sequence.
// 
  else if ( q[0] == m )
  {
    more = false;

    for ( i = 0; i < k - 1; i++ )
    {
      q[i] = 0;
    }
    q[k-1] = m;

    leftmost = k + 1;
  }
  else if ( leftmost < k + 1 )
  {
    leftmost = leftmost - 1;
    q[k-1] = q[leftmost-1] - 1;
    q[leftmost-1] = 0;
    q[leftmost-2] = q[leftmost-2] + 1;
    if ( q[k-1] != 0 )
    {
      leftmost = k + 1;
    }
  }
  else
  {
    if ( q[k-1] == 1 )
    {
      leftmost = k;
    }
    q[k-1] = q[k-1] - 1;
    q[k-2] = q[k-2] + 1;
  }
  return;
}
//****************************************************************************80

int edge_check ( int n_node, int n_edge, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    EDGE_CHECK checks a graph stored by edges.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N_NODE, the number of nodes in the graph.
//    N_NODE must be positive.
// 
//    Input, int N_EDGE, the number of edges in the graph.
//    N_EDGE must be positive.
// 
//    Input, int T(2,N_EDGE), describes the edges of the tree
//    as pairs of nodes.
// 
//    Output, int EDGE_CHECK, error flag.
//    -1, N_NODE is not positive.
//    -2, N_EDGE is not positive.
//    0, no error.
//    I, edge T(1,I), T(2,I) is illegal.
// 
{
  int i;
  int ierror;
  int j;
  int j2;

  ierror = 0;

  if ( n_node < 1 )
  {
    cerr << "\n";
    cerr << "EDGE_CHECK - Fatal error!\n";
    cerr << "  N < 1\n";
    exit ( 1 );
  }

  if ( n_edge < 1 )
  {
    ierror = -2;
    return ierror;
  }
// 
//  Every edge must join two legal nodes.
// 
  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n_edge; j++ )
    {
      if ( t[i+j*2] < 1 || n_node < t[i+j*2] )
      {
        ierror = j + 1;
        return ierror;
      }
    }
  }
// 
//  Every edge must join distinct nodes.
// 
  for ( j = 0; j < n_edge; j++ )
  {
    if ( t[0+j*2] == t[1+j*2] )
    {
      ierror = j + 1;
      return ierror;
    }
  }
// 
//  Every edge must be distinct.
// 
  for ( j = 0; j < n_edge - 1; j++ )
  {
    for ( j2 = j + 1; j2 < n_edge; j2++ )
    {
      if ( t[0+j*2] == t[0+j2*2] && t[1+j*2] == t[1+j2*2] )
      {
        ierror = j2 + 1;
        return ierror;
      }
      else if ( t[0+j*2] == t[1+j2*2] && t[1+j*2] == t[0+j2*2] )
      {
        ierror = j2 + 1;
        return ierror;
      }
    }
  }
  return ierror;
}
//****************************************************************************80

int *edge_degree ( int n_node, int n_edge, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    EDGE_DEGREE returns the degree of the nodes of a graph stored by edges.
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
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
// 
//  Parameters:
// 
//    Input, int N_NODE, the number of nodes in the graph.
//    N_NODE must be positive.
// 
//    Input, int N_EDGE, the number of edges in the graph.
//    N_EDGE must be positive.
// 
//    Input, int T[2*N_EDGE], describes the edges of the tree
//    as pairs of nodes.
// 
//    Output, int EDGE_DEGREE[N_NODE], the degree of each node.
// 
{
  int *d;
  int i;
  int ierror;
  int j;
// 
//  Check.
// 
  ierror = edge_check ( n_node, n_edge, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "EDGE_DEGREE - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Compute the degree of each node.
//
  d = new int[n_node];

  for ( i = 0; i < n_node; i++ )
  {
    d[i] = 0;
  }
  for ( j = 0; j < n_edge; j++ )
  {
    d[t[0+j*2]-1] = d[t[0+j*2]-1] + 1;
    d[t[1+j*2]-1] = d[t[1+j*2]-1] + 1;
  }

  return d;
}
//****************************************************************************80

int edge_enum ( int n_node )

//****************************************************************************80
// 
//  Purpose:
//
//    EDGE_ENUM enumerates the maximum number of edges in a graph on N_NODE nodes.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N_NODE, the number of nodes in the graph.
//    N_NODE must be positive.
// 
//    Output, int EDGE_ENUM, the maximum number of edges in a graph
//    on N_NODE nodes.
// 
{
  int value;

  value = ( n_node * ( n_node - 1 ) ) / 2;

  return value;
}
//****************************************************************************80

int fall ( int x, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    FALL computes the falling factorial function [X]_N.
// 
//  Discussion:
// 
//    The number of "injections" or 1-to-1 mappings from
//    a set of N elements to a set of M elements is [M]_N.
// 
//    The number of permutations of N objects out of M is [M}_N.
// 
//    The Stirling numbers of the first kind can be used
//    to convert a falling factorial into a polynomial, as follows:
// 
//      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
// 
//  Formula:
// 
//    [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int X, the argument of the falling factorial
//    function.
// 
//    Input, int N, the order of the falling factorial function.
//    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
//    negative, a "rising" factorial will be computed.
// 
//    Output, int FALL, the falling factorial function.
// 
{
  int arg;
  int i;
  int value;

  value = 1;

  arg = x;

  if ( 0 < n )
  {
    for ( i = 1; i <= n; i++ )
    {
      value = value * arg;
      arg = arg - 1;
    }
  }
  else if ( n < 0 )
  {
    for ( i = n; i <= -1; i++ )
    {
      value = value * arg;
      arg = arg + 1;
    }
  }
  return value;
}
//****************************************************************************80

int gray_code_check ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    GRAY_CODE_CHECK checks a Gray code element.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of digits in each element.
//    N must be positive.
// 
//    Input, int T[N], an element of the Gray code.
//    Each entry T(I) is either 0 or 1.
// 
//    Output, int GRAY_CODE_CHECK, error flag.
//    0, no error, T represents a Gray code element.
//    -1, N is not positive.
//    I, error, T(I) is an illegal value for a Gray code element.
// 
{
  int i;
  int ierror;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] != 0 && t[i] != 1 )
    {
      ierror = i + 1;
      return ierror;
    }
  }

  return ierror;
}
//****************************************************************************80

int gray_code_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    GRAY_CODE_ENUM enumerates the Gray codes on N digits.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of digits in each element.
//    N must be nonnegative.
// 
//    Output, int GRAY_CODE_ENUM, the number of distinct elements.
// 
{
  int value;

  value = i4_power ( 2, n );

  return value;
}
//****************************************************************************80

int gray_code_rank ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    GRAY_CODE_RANK computes the rank of a Gray code element.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the number of digits in each element.
//    N must be positive.
// 
//    Input, int T[N], an element of the Gray code.
//    Each entry is either 0 or 1.
// 
//    Output, int GRAY_CODE_RANK, the rank of the element.
// 
{
  int b;
  int i;
  int ierror;
  int rank;
// 
//  Check.
// 
  ierror = gray_code_check ( n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "GRAY_CODE_RANK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  rank = 0;
  b = 0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( t[n-i-1] != 0 )
    {
      b = 1 - b;
    }
    if ( b == 1 )
    {
      rank = rank + i4_power ( 2, i );
    }
  }
  return rank;
}
//****************************************************************************80

void gray_code_successor ( int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    GRAY_CODE_SUCCESSOR computes the binary reflected Gray code successor.
// 
//  Example:
// 
//    000, 001, 011, 010, 110, 111, 101, 100,
//    after which the sequence repeats.
// 
//  Discussion:
// 
//    In the original code, the successor of the element that has an
//    initial 1 followed by N-1 zeroes is undefined.  In this version,
//    the successor is the element with N zeroes.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of digits in each element.
//    N must be positive.
// 
//    Input/output, int T[N].
//    On input, T contains an element of the Gray code, that is,
//    each entry T(I) is either 0 or 1.
//    On output, T contains the successor to the input value; this
//    is an element of the Gray code, which differs from the input
//    value in a single position.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
  int weight;
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
// 
//  Check.
// 
  ierror = gray_code_check ( n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "GRAY_CODE_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  weight = i4vec_sum ( n, t );

  if ( ( weight % 2 ) == 0 )
  {
    if ( t[n-1] == 0 )
    {
      t[n-1] = 1;
    }
    else
    {
      t[n-1] = 0;
    }
    rank = rank + 1;
    return;
  }
  else
  {
    for ( i = n - 1; 1 <= i; i-- )
    {
      if ( t[i] == 1 )
      {
        if ( t[i-1] == 0 )
        {
          t[i-1] = 1;
        }
        else
        {
          t[i-1] = 0;
        }
        rank = rank + 1;
        return;
      }
    }
// 
//  The final element was input.
//  Return the first element.
// 
    for ( i = 0; i < n; i++ )
    {
      t[i] = 0;
    }
    rank = 0;
  }
  return;
}
//****************************************************************************80

int *gray_code_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    GRAY_CODE_UNRANK computes the Gray code element of given rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the element.
//    0 <= RANK <= 2^N.
// 
//    Input, int N, the number of digits in each element.
//    N must be positive.
// 
//    Output, int GRAY_CODE_UNRANK[N], the element of the Gray code which has
//    the given rank.
// 
{
  int b;
  int bprime;
  int i;
  int ngray;
  int rank_copy;
  int *t;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "GRAY_CODE_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  ngray = gray_code_enum ( n );

  if ( rank < 0 || ngray < rank )
  {
    cout << "\n";
    cout << "GRAY_CODE_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  t = new int[n];

  rank_copy = rank;
  for ( i = 0; i < n; i++ )
  {
    t[i] = 0;
  }
  bprime = 0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b = rank_copy / i4_power ( 2, i );

    if ( b != bprime )
    {
      t[n-i-1] = 1;
    }
    bprime = b;
    rank_copy = rank_copy - b * i4_power ( 2, i );
  }
  return t;
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

int i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL computes the factorial of N.
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
//    26 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    If N is less than 1, the function value is returned as 1.
//    0 <= N <= 13 is required.
//
//    Output, int I4_FACTORIAL, the factorial of N.
//
{
  int i;
  int value;

  value = 1;

  if ( 13 < n )
  {
    cerr << "I4_FACTORIAL - Fatal error!\n";
    cerr << "  I4_FACTORIAL(N) cannot be computed as an integer\n";
    cerr << "  for 13 < N.\n";
    cerr << "  Input value N = " << n << "\n";
    exit ( 1 );
  }

  for ( i = 1; i <= n; i++ )
  {
    value = value * i;
  }

  return value;
}
//****************************************************************************80

void i4_factorial_values ( int &n_data, int &n, int &fn )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL_VALUES returns values of the factorial function.
//
//  Discussion:
//
//    0! = 1
//    I! = Product ( 1 <= J <= I ) I
//
//    In Mathematica, the function can be evaluated by:
//
//      n!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2011
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
//    Output, int &N, the argument of the function.
//
//    Output, int &FN, the value of the function.
//
{
# define N_MAX 13

  static int fn_vec[N_MAX] = {
            1,
            1,
            2,
            6,
           24,
          120,
          720,
         5040,
        40320,
       362880,
      3628800,
     39916800,
    479001600 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
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

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
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
//    12 November 2006
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
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

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
//  Print the columns of the matrix, in strips of INCX.
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
    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
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
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_backtrack ( int n, int maxstack, int stack[], int x[], int *indx, 
  int *k, int *nstack, int ncan[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
//
//  Discussion:
//
//    The routine tries to construct an integer vector one index at a time,
//    using possible candidates as supplied by the user.
//
//    At any time, the partially constructed vector may be discovered to be
//    unsatisfactory, but the routine records information about where the
//    last arbitrary choice was made, so that the search can be
//    carried out efficiently, rather than starting out all over again.
//
//    First, call the routine with INDX = 0 so it can initialize itself.
//
//    Now, on each return from the routine, if INDX is:
//      1, you've just been handed a complete candidate vector;
//         Admire it, analyze it, do what you like.
//      2, please determine suitable candidates for position X(K).
//         Return the number of candidates in NCAN(K), adding each
//         candidate to the end of STACK, and increasing NSTACK.
//      3, you're done.  Stop calling the routine;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2004
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
//    Input, int N, the number of positions to be filled in the vector.
//
//    Input, int MAXSTACK, the maximum length of the stack.
//
//    Input, int STACK[MAXSTACK], a list of all current candidates for
//    all positions 1 through K.
//
//    Input/output, int X[N], the partial or complete candidate vector.
//
//    Input/output, int *INDX, a communication flag.
//    On input,
//      0 to start a search.
//    On output:
//      1, a complete output vector has been determined and returned in X(1:N);
//      2, candidates are needed for position X(K);
//      3, no more possible vectors exist.
//
//    Input/output, int *K, if INDX=2, the current vector index being considered.
//
//    Input/output, int *NSTACK, the current length of the stack.
//
//    Input/output, int NCAN[N], lists the current number of candidates for
//    positions 1 through K.
//
{
//
//  If this is the first call, request a candidate for position 1.
//
  if ( *indx == 0 )
  {
    *k = 1;
    *nstack = 0;
    *indx = 2;
    return;
  }
//
//  Examine the stack.
//
  for ( ; ; )
  {
//
//  If there are candidates for position K, take the first available
//  one off the stack, and increment K.
//
//  This may cause K to reach the desired value of N, in which case
//  we need to signal the user that a complete set of candidates
//  is being returned.
//
    if ( 0 < ncan[(*k)-1] )
    {
      x[(*k)-1] = stack[(*nstack)-1];
      *nstack = *nstack - 1;

      ncan[(*k)-1] = ncan[(*k)-1] - 1;

      if ( *k != n )
      {
        *k = *k + 1;
        *indx = 2;
      }
      else
      {
        *indx = 1;
      }
      break;
    }
//
//  If there are no candidates for position K, then decrement K.
//  If K is still positive, repeat the examination of the stack.
//
    else
    {
      *k = *k - 1;

      if ( *k <= 0 )
      {
        *indx = 3;
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

void i4vec_indicator ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int A[N], the initialized array.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return;
}
//****************************************************************************80

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
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
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MAX, the value of the maximum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
//****************************************************************************80

int *i4vec_part1 ( int n, int npart )

//****************************************************************************80
// 
//  Purpose:
//
//    I4VEC_PART1 partitions an integer N into NPART parts.
// 
//  Example:
// 
//    Input:
// 
//      N = 17, NPART = 5
// 
//    Output:
// 
//      X = ( 13, 1, 1, 1, 1 ).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the integer to be partitioned.  N
//    may be positive, zero, or negative.
// 
//    Input, int NPART, the number of entries in the array.
//    1 <= NPART <= N.
// 
//    Output, int I4VEC_PART1[NPART], the partition of N.  The entries of
//    X add up to N.  X(1) = N + 1 - NPART, and all other entries
//    are equal to 1.
// 
{
  int i;
  int *x;

  if ( npart < 1 || n < npart )
  {
    cout << "\n";
    cout << "I4VEC_PART1 - Fatal error!\n";
    cout << "  The input value of NPART is illegal.\n";
    exit ( 1 );
  }

  x = new int[npart];

  x[0] = n + 1 - npart;
  for ( i = 1; i < npart; i++ )
  {
    x[i] = 1;
  }

  return x;
}
//****************************************************************************80

void i4vec_part2 ( int n, int npart, int x[] )

//****************************************************************************80
// 
//  Purpose:
//
//    I4VEC_PART2 partitions an integer N into NPART nearly equal parts.
// 
//  Discussion:
//
//    Thanks to John Nitao for pointing out a typographical error in 
//    of the form "x[j=1]" for "x[j-1]", 14 April 2013.
//
//  Example:
// 
//    Input:
// 
//      N = 17, NPART = 5
// 
//    Output:
// 
//      X = ( 4, 4, 3, 3, 3 ).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    14 April 2013
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the integer to be partitioned.  N
//    may be positive, zero, or negative.
// 
//    Input, int NPART, the number of entries in the array.
//    1 <= NPART
// 
//    Output, int X[NPART], the partition of N.  The entries of
//    X add up to N.  The entries of X are either all equal, or
//    differ by at most 1.  The entries of X all have the same sign
//    as N, and the "largest" entries occur first.
// 
{
  int i;
  int j;

  if ( npart < 1 )
  {
    cout << "\n";
    cout << "I4VEC_PART2 - Fatal error!\n";
    cout << "  The input value of NPART is illegal.\n";
    exit ( 1 );
  }

  for ( i = 0; i < npart; i++ )
  {
    x[i] = 0;
  }

  if ( 0 < n )
  {
    j = 1;
    for ( i = 1; i <= n; i++ )
    {
      x[j-1] = x[j-1] + 1;
      j = j + 1;
      if ( npart < j )
      {
        j = 1;
      }
    }
  }
  else if ( n < 0 )
  {
    j = 1;
    for ( i = n; i <= -1; i++ )
    {
      x[j-1] = x[j-1] - 1;
      j = j + 1;
      if ( npart < j )
      {
        j = 1;
      }
    }
  }

  return;
}
//****************************************************************************80

int *i4vec_part2_new ( int n, int npart )

//****************************************************************************80
// 
//  Purpose:
//
//    I4VEC_PART2_NEW partitions an integer N into NPART nearly equal parts.
// 
//  Example:
// 
//    Input:
// 
//      N = 17, NPART = 5
// 
//    Output:
// 
//      X = ( 4, 4, 3, 3, 3 ).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the integer to be partitioned.  N
//    may be positive, zero, or negative.
// 
//    Input, int NPART, the number of entries in the array.
//    1 <= NPART
// 
//    Output, int I4VEC_PART2[NPART], the partition of N.  The entries of
//    X add up to N.  The entries of X are either all equal, or
//    differ by at most 1.  The entries of X all have the same sign
//    as N, and the "largest" entries occur first.
// 
{
  int *x;

  x = new int[npart];

  i4vec_part2 ( n, npart, x );

  return x;
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

void i4vec_reverse ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N], the array to be reversed.
//
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}
//****************************************************************************80

int i4vec_search_binary_a ( int n, int a[], int b )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    Binary search is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.9,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 26.
//
//  Parameters:
//
//    Input, int N, the number of elements in the vector.
//
//    Input, int A[N], the array to be searched.  A must
//    be sorted in ascending order.
//
//    Input, int B, the value to be searched for.
//
//    Output, int I4VEC_SEARCH_BINARY_A, the result of the search.
//    -1, B does not occur in A.
//    I, A[I] = B.
//
{
  int high;
  int index;
  int low;
  int mid;
//
//  Check.
//
  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_SEARCH_BINARY_A - Fatal error!\n";
    cerr << "  The array dimension N is less than 1.\n";
    exit ( 1 );
  }

  index = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == b )
    {
      index = mid;
      break;
    }
    else if ( a[mid-1] < b )
    {
      low = mid + 1;
    }
    else if ( b < a[mid-1] )
    {
      high = mid - 1;
    }
  }
  return index;
}
//****************************************************************************80

int i4vec_search_binary_d ( int n, int a[], int b )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC for a value.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    Binary search is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.9,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 26.
//
//  Parameters:
//
//    Input, int N, the number of elements in the vector.
//
//    Input, int A[N], the array to be searched.  A must
//    be sorted in descending order.
//
//    Input, int B, the value to be searched for.
//
//    Output, int I4VEC_SEARCH_BINARY_D, the result of the search.
//    -1, B does not occur in A.
//    I, A[I] = B.
//
{
  int high;
  int index;
  int low;
  int mid;
//
//  Check.
//
  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_SEARCH_BINARY_D - Fatal error!\n";
    cerr << "  The array dimension N is less than 1.\n";
    exit ( 1 );
  }

  index = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == b )
    {
      index = mid;
      break;
    }
    else if ( b < a[mid-1] )
    {
      low = mid + 1;
    }
    else if ( a[mid-1] < b )
    {
      high = mid - 1;
    }
  }
  return index;
}
//****************************************************************************80

void i4vec_sort_insert_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_INSERT_A uses an ascending insertion sort on an I4VEC.
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
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.1,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 11.
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//    N must be positive.
//
//    Input/output, int A[N].
//    On input, A contains data to be sorted.
//    On output, the entries of A have been sorted in ascending order.
//
{
  int i;
  int j;
  int x;

  for ( i = 1; i < n; i++ )
  {
    x = a[i];

    j = i;

    while ( 1 <= j && x < a[j-1] )
    {
      a[j] = a[j-1];
      j = j - 1;
    }

    a[j] = x;
  }

  return;
}
//****************************************************************************80

void i4vec_sort_insert_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_INSERT_D uses a descending insertion sort on an I4VEC.
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
//    13 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Kreher, Douglas Simpson,
//    Algorithm 1.1,
//    Combinatorial Algorithms,
//    CRC Press, 1998, page 11.
//
//  Parameters:
//
//    Input, int N, the number of items in the vector.
//    N must be positive.
//
//    Input/output, int A[N].
//    On input, A contains data to be sorted.
//    On output, the entries of A have been sorted in ascending order.
//
{
  int i;
  int j;
  int x;

  for ( i = 1; i < n; i++ )
  {
    x = a[i];
    j = i;

    while ( 1 <= j && a[j-1] < x )
    {
      a[j] = a[j-1];
      j = j - 1;
    }
    a[j] = x;
  }

  return;
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

void i4vec_transpose_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
//    TITLE = "My vector:  "
//
//    My vector:      1    2    3    4    5
//                    6    7    8    9   10
//                   11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2004
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
  int ihi;
  int ilo;
  int title_len;

  title_len = title.length ( );

  for ( ilo = 1; ilo <= n; ilo = ilo + 10 )
  {
    ihi = i4_min ( ilo + 10 - 1, n );
    if ( ilo == 1 )
    {
      cout << title;
    }
    else
    {
      for ( i = 1; i <= title_len; i++ )
      {
        cout << " ";
      }
    }
    for ( i = ilo; i <= ihi; i++ )
    {
      cout << "  " << setw(5) << a[i-1];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void knapsack_01 ( int n, double mass_limit, double p[], double w[], double x[], 
  double &mass, double &profit )

//****************************************************************************80
// 
//  Purpose:
//
//    KNAPSACK_01 solves the 0/1 knapsack problem.
// 
//  Discussion:
// 
//    The 0/1 knapsack problem is as follows:
// 
//      Given:
//        a set of N objects,
//        a profit P(I) and weight W(I) associated with each object,
//        and a weight limit MASS_LIMIT,
//      Determine:
//        a set of choices X(I) which are 0 or 1, that maximizes the profit
//          P = Sum ( 1 <= I <= N ) P(I) * X(I)
//        subject to the constraint
//          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
// 
//    This routine assumes that the objects have already been sorted
//    in order of decreasing "profit density", P(I)/W(I).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int N, the number of objects.
// 
//    Input, double MASS_LIMIT, the weight limit of the
//    chosen objects.
// 
//    Input/output, double P[N], the "profit" or value of each object.
//    P is assumed to be nonnegative.
// 
//    Input/output, double W[N], the "weight" or cost of each object.
//    W is assumed to be  nonnegative.
// 
//    Output, double X[N], the choice function for the objects.
//    0, the object was not taken.
//    1, the object was taken.
// 
//    Output, double &MASS, the total mass of the objects taken.
// 
//    Output, double &PROFIT, the total profit of the objects taken.
// 
{
  int i;
  int indx;
  int k;
  double mass_1;
  double mass_2;
  double mass_best;
  double mass_remaining;
  int maxstack = 100;
  int *ncan;
  int nstack;
  double profit_1;
  double profit_2;
  double profit_best;
  double *stack;
  double *x_best;

  ncan = new int[n];
  stack = new double[maxstack];
  x_best = new double[n];

  nstack = 0;
// 
//  Initialize the "best so far" data.
// 
  for ( i = 0; i < n; i++ )
  {
    x_best[i] = 0.0;
  }
  profit_best = 0.0;
  mass_best = 0;
// 
//  Begin the backtracking solution.
// 
  indx = 0;

  for ( ; ; )
  {
    r8vec_backtrack ( n, maxstack, stack, x, &indx, &k, &nstack, ncan );
// 
//  Got a new candidate.  Compare it to the best so far.
// 
    if ( indx == 1 )
    {
      profit = r8vec_dot_product ( n, p, x );
      mass = r8vec_dot_product ( n, w, x );

      if ( profit_best < profit || ( profit == profit_best && mass < mass_best ) )
      {
        profit_best = profit;
        mass_best = mass;
        for ( i = 0; i < n; i++ )
        {
          x_best[i] = x[i];
        }
      }
    }
// 
//  Need candidates for X(K).
// 
//  X(K) = 1 is possible if:
// 
//    * adding W(K) to our mass doesn''t put us over our mass limit;
//    * and adding P(K) to our current profit, and taking the best we
//      could get using rational X for the remainder would put us over
//      our current best.
// 
//  X(K) = 0 is always possible.
// 
    else if ( indx == 2 )
    {
      ncan[k-1] = 0;

      mass_1 = w[k-1];
      for ( i = 0; i < k - 1; i++ )
      {
        mass_1 = mass_1 + w[i] * x[i];
      }

      if ( mass_1 <= mass_limit )
      {
        mass_remaining = mass_limit - mass_1;

        profit_1 = p[k-1];
        for ( i = 0; i < k - 1; i++ )
        {
          profit_1 = profit_1 + p[i] * x[i];
        }

        if ( k < n )
        {
          knapsack_rational ( n - k, mass_remaining, p+k, w+k, 
            x+k, mass_2, profit_2 );
        }
        else
        {
          profit_2 = 0.0;
        }

        if ( profit_best < profit_1 + profit_2 )
        {
          if ( maxstack <= nstack )
          {
            cout << "\n";
            cout << "KNAPSACK_01 - Fatal error!\n";
            cout << "  Exceeded stack space.\n";
            return;
          }
          ncan[k-1] = ncan[k-1] + 1;
          nstack = nstack + 1;
          stack[nstack-1] = 1.0;
        }
      }

      if ( maxstack <= nstack )
      {
        cout << "\n";
        cout << "KNAPSACK_01 - Fatal error!\n";
        cout << "  Exceeded stack space.\n";
        return;
      }

      ncan[k-1] = ncan[k-1] + 1;
      nstack = nstack + 1;
      stack[nstack-1] = 0.0;
    }
// 
//  Done.  Return the best solution.
// 
    else
    {
      profit = profit_best;
      mass = mass_best;
      for ( i = 0; i < n; i++ )
      {
        x[i] = x_best[i];
      }
      break;
    }
  }

  delete [] ncan;
  delete [] stack;
  delete [] x_best;

  return;
}
//****************************************************************************80

void knapsack_rational ( int n, double mass_limit, double p[], double w[],
  double x[], double &mass, double &profit )

//****************************************************************************80
// 
//  Purpose:
//
//    KNAPSACK_RATIONAL solves the rational knapsack problem.
// 
//  Discussion:
// 
//    The rational knapsack problem is a generalization of the 0/1 knapsack
//    problem.  It is mainly used to derive a bounding function for the
//    0/1 knapsack problem.
// 
//    The 0/1 knapsack problem is as follows:
// 
//      Given:
//        a set of N objects,
//        a profit P(I) and weight W(I) associated with each object,
//        and a weight limit MASS_LIMIT,
//      Determine:
//        a set of choices X(I) which are 0 or 1, that maximizes the profit
//          P = Sum ( 1 <= I <= N ) P(I) * X(I)
//        subject to the constraint
//          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
// 
//    By contrast, the rational knapsack problem allows the values X(I)
//    to be any value between 0 and 1.  A solution for the rational knapsack
//    problem is known.  Arrange the objects in order of their "profit density"
//    ratios P(I)/W(I), and then take in order as many of these as you can.
//    If you still have "room" in the weight constraint, then you should
//    take the maximal fraction of the very next object, which will complete
//    your weight limit, and maximize your profit.
// 
//    If should be obvious that, given the same data, a solution for
//    the rational knapsack problem will always have a profit that is
//    at least as high as for the 0/1 problem.  Since the rational knapsack
//    maximum profit is easily computed, this makes it a useful bounding
//    function.
// 
//    Note that this routine assumes that the objects have already been
//    arranged in order of the "profit density".
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of objects.
// 
//    Input, double MASS_LIMIT, the weight limit of the
//    chosen objects.
// 
//    Input, double P[N], the "profit" or value of each object.
//    The entries of P are assumed to be nonnegative.
// 
//    Input, double W[N], the "weight" or cost of each object.
//    The entries of W are assumed to be nonnegative.
// 
//    Output, double X[N], the choice function for the objects.
//    0.0, the object was not taken.
//    1.0, the object was taken.
//    R, where 0 < R < 1, a fractional amount of the object was taken.
// 
//    Output, double &MASS, the total mass of the objects taken.
// 
//    Output, double &PROFIT, the total profit of the objects taken.
// 
{
  int i;

  mass = 0.0;
  profit = 0.0;

  for ( i = 0; i < n; i++ )
  {
    if ( mass_limit <= mass )
    {
      x[i] = 0.0;
    }
    else if ( mass + w[i] <= mass_limit )
    {
      x[i] = 1.0;
      mass = mass + w[i];
      profit = profit + p[i];
    }
    else
    {
      x[i] = ( mass_limit - mass ) / w[i];
      mass = mass_limit;
      profit = profit + p[i] * x[i];
    }
  }
  return;
}
//****************************************************************************80

void knapsack_reorder ( int n, double p[], double w[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KNAPSACK_REORDER reorders the knapsack data by "profit density".
// 
//  Discussion:
// 
//    This routine must be called to rearrange the data before calling
//    routines that handle a knapsack problem.
// 
//    The "profit density" for object I is defined as P(I)/W(I).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of objects.
// 
//    Input/output, double P[N], the "profit" or value of each object.
// 
//    Input/output, double W[N], the "weight" or cost of each object.
// 
{
  int i;
  int j;
  double t;
// 
//  Rearrange the objects in order of "profit density".
//
  for ( i = 0; i < n; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      if ( p[i] * w[j] < p[j] * w[i] )
      {
        t    = p[i];
        p[i] = p[j];
        p[j] = t;

        t    = w[i];
        w[i] = w[j];
        w[j] = t;
      }
    }
  }
  return;
}
//****************************************************************************80

int ksubset_colex_check ( int k, int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_COLEX_CHECK checks a K subset in colex form.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have. 1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input, int T[K], describes a K subset.  T(I) is the I-th
//    element of the K subset.  The elements must be listed in
//    DESCENDING order.
// 
//    Output, int KSUBSET_COLEX_CHECK, error flag.
//    0, no error.
//    -1, N is not positive.
//    -2, K is not positive.
//    I, entry I is illegal.
// 
{
  int i;
  int ierror;
  int tmax;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  if ( k < 1 || n < k )
  {
    ierror = -2;
    return ierror;
  }

  tmax = n + 1;

  for ( i = 0; i < k; i++ )
  {
    if ( t[i] <= 0 || tmax <= t[i] )
    {
      ierror = i;
      return ierror;
    }
    tmax = t[i];
  }
  return ierror;
}
//****************************************************************************80

int ksubset_colex_rank ( int k, int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_COLEX_RANK computes the colex rank of a K subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input, int T[K[, describes a K subset.  T(I) is the I-th
//    element of the K subset.  The elements must be listed in DESCENDING order.
// 
//    Output, int KSUBSET_COLEX_RANK, the rank of the subset.
// 
{
  int i;
  int ierror;
  int rank;
// 
//  Check.
// 
  ierror = ksubset_colex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_COLEX_CHECK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  rank = 0;

  for ( i = 0; i < k; i++ )
  {
    rank = rank + i4_choose ( t[i] - 1, k - i );
  }

  return rank;
}
//****************************************************************************80

void ksubset_colex_successor ( int k, int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_COLEX_SUCCESSOR computes the K subset colex successor.
// 
//  Discussion:
// 
//    In the original code, there is a last element with no successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input/output, int T[K], describes a K subset.  T(I) is the
//    I-th element.  The elements must be listed in DESCENDING order.
//    On input, T describes a K subset.
//    On output, T describes the next K subset in the ordering.
//    If the input T was the last in the ordering, then the output T
//    will be the first.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    for ( i = 1; i <= k; i++ )
    {
      t[i-1] = k + 1 - i;
    }
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = ksubset_colex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_COLEX_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  for ( i = k - 1; 1 <= i; i-- )
  {
    if ( t[k-i] + 1 < t[k-i-1] )
    {
      t[k-i] = t[k-i] + 1;
      rank = rank + 1;
      return;
    }
  }

  if ( t[0] < n )
  {
    t[0] = t[0] + 1;
    for ( i = 1; i <= k - 1; i++ )
    {
      t[k-i] = i;
    }
    rank = rank + 1;
    return;
  }
// 
//  The last K subset was input.
//  Return the first one.
// 
  for ( i = 1; i <= k; i++ )
  {
    t[i-1] = k + 1 - i;
  }

  rank = 0;

  return;
}
//****************************************************************************80

int *ksubset_colex_unrank ( int rank, int k, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_COLEX_UNRANK computes the K subset of given colex rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the K subset.
// 
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Output, int KSUBSET_COLEX_UNRANK[K], describes the K subset of the given
//    rank.  T(I) is the I-th element.  The elements must be listed in
//    DESCENDING order.
// 
{
  int i;
  int nksub;
  int rank_copy;
  int *t;
  int x;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "KSUBSET_COLEX_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  if ( k < 1 || n < k )
  {
    cout << "\n";
    cout << "KSUBSET_COLEX_UNRANK - Fatal error!\n";
    cout << "  Input K is illegal.\n";
    exit ( 1 );
  }

  nksub = ksubset_enum ( k, n );

  if ( rank < 0 || nksub < rank )
  {
    cout << "\n";
    cout << "KSUBSET_COLEX_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }
// 
  rank_copy = rank;

  x = n;

  t = new int[k];

  for ( i = 1; i <= k; i++ )
  {
    while ( rank_copy < i4_choose ( x, k + 1 - i ) )
    {
      x = x - 1;
    }

    t[i-1] = x + 1;
    rank_copy = rank_copy - i4_choose ( x, k + 1 - i );
  }

  return t;
}
//****************************************************************************80

int ksubset_enum ( int k, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_ENUM enumerates the K element subsets of an N set.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int K, the number of elements each K subset must
//    have. 0 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    0 <= N.
// 
//    Output, int KSUBSET_ENUM, the number of distinct elements.
// 
{
  int value;

  value = i4_choose ( n, k );

  return value;
}
//****************************************************************************80

int ksubset_lex_check ( int k, int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_LEX_CHECK checks a K subset in lex form.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have. 1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input, int T[K], describes a K subset.  T(I) is the I-th
//    element of the K subset.  The elements must be listed in
//    DESCENDING order.
// 
//    Output, int IERROR, error flag.
//    0, no error.
//    -1, N is illegal.
//    -2, K is illegal.
//    I+1, entry I is illegal.
// 
{
  int i;
  int ierror;
  int tmin;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  if ( k < 1 || n < k )
  {
    ierror = -2;
    return ierror;
  }

  tmin = 0;

  for ( i = 0; i < k; i++ )
  {
    if ( t[i] <= tmin || n < t[i] )
    {
      ierror = i + 1;
      return ierror;
    }
    tmin = t[i];
  }
  return ierror;
}
//****************************************************************************80

int ksubset_lex_rank ( int k, int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_LEX_RANK computes the lexicographic rank of a K subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input, int T[K], describes a K subset.  T(I) is the I-th
//    element.  The elements must be listed in ascending order.
// 
//    Output, int KSUBSET_LEX_RANK, the rank of the K subset.
// 
{
  int i;
  int ierror;
  int j;
  int rank;
  int tim1;
// 
//  Check.
// 
  ierror = ksubset_lex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_LEX_RANK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    exit ( 1 );
  }

  rank = 0;

  for ( i = 1; i <= k; i++ )
  {
    if ( i == 1 )
    {
      tim1 = 0;
    }
    else
    {
      tim1 = t[i-2];
    }

    if ( tim1 + 1 <= t[i-1] - 1 )
    {
      for ( j = tim1 + 1; j <= t[i-1] - 1; j++ )
      {
        rank = rank + i4_choose ( n - j, k - i );
      }
    }
  }

  return rank;
}
//****************************************************************************80

void ksubset_lex_successor ( int k, int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_LEX_SUCCESSOR computes the K subset lexicographic successor.
// 
//  Discussion:
// 
//    In the original code, there is a last element with no successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have. 1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input/output, int T[K], describes a K subset.  T(I) is
//    the I-th element.  The elements must be listed in ascending order.
//    On input, T describes a K subset.
//    On output, T describes the next K subset in the ordering.
//    If the input T was the last in the ordering, then the output T
//    will be the first.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
  int isave;
  int j;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    i4vec_indicator ( k, t );
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = ksubset_lex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_LEX_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  isave = 0;

  for ( i = k; 1 <= i; i-- )
  {
    if ( t[i-1] != n - k + i )
    {
      isave = i;
      break;
    }
  }
// 
//  The last K subset was input.
//  Return the first one.
// 
  if ( isave == 0 )
  {
    i4vec_indicator ( k, t );
    rank = 0;
  }
  else
  {
    for ( j = k; isave <= j; j-- )
    {
      t[j-1] = t[isave-1] + 1 + j - isave;
    }
    rank = rank + 1;
  }

  return;
}
//****************************************************************************80

int *ksubset_lex_unrank ( int rank, int k, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_LEX_UNRANK computes the K subset of given lexicographic rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the K subset.
// 
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Output, int KSUBSET_LEX_RANK[K], describes the K subset of the given
//    rank.  T(I) is the I-th element.  The elements must be listed in
//    ascending order.
// 
{
  int i;
  int nksub;
  int rank_copy;
  int *t;
  int x;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "KSUBSET_LEX_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  if ( k < 1 || n < k )
  {
    cout << "\n";
    cout << "KSUBSET_LEX_UNRANK - Fatal error!\n";
    cout << "  Input K is illegal.\n";
    exit ( 1 );
  }

  nksub = ksubset_enum ( k, n );

  if ( rank < 0 || nksub < rank )
  {
    cout << "\n";
    cout << "KSUBSET_LEX_UNRANK - Fatal error!\n";
    cout << "  Input rank is illegal.\n";
    exit ( 1 );
  }

  t = new int[k];

  rank_copy = rank;

  x = 1;

  for ( i = 1; i <= k; i++ )
  {
    while ( i4_choose ( n - x, k - i ) <= rank_copy )
    {
      rank_copy = rank_copy - i4_choose ( n - x, k - i );
      x = x + 1;
    }

    t[i-1] = x;
    x = x + 1;
  }

  return t;
}
//****************************************************************************80

int ksubset_revdoor_rank ( int k, int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_REVDOOR_RANK computes the revolving door rank of a K subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input, int T[K], describes a K subset.  T(I) is the I-th
//    element.  The elements must be listed in ascending order.
// 
//    Output, int KSUBSET_REVDOOR_RANK, the rank of the K subset.
// 
{
  int i;
  int ierror;
  int rank;
  int s;
// 
//  Check.
// 
  ierror = ksubset_lex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_REVDOOR_RANK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  if ( ( k % 2 ) == 0 )
  {
    rank = 0;
  }
  else
  {
    rank = - 1;
  }

  s = 1;

  for ( i = k; 1 <= i; i-- )
  {
    rank = rank + s * i4_choose ( t[i-1], i );
    s = - s;
  }

  return rank;
}
//****************************************************************************80

void ksubset_revdoor_successor ( int k, int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_REVDOOR_SUCCESSOR computes the K subset revolving door successor.
// 
//  Discussion:
// 
//    After numerous attempts to implement the algorithm published in
//    Kreher and Stinson, the Nijenhuis and Wilf version was implemented
//    instead.  The K and S algorithm is supposedly based on the N and W one.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
// 
//  Parameters:
// 
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Input/output, int T[K], describes a K subset.  T(I) is the
//    I-th element.  The elements must be listed in ascending order.
//    On input, T describes a K subset.
//    On output, T describes the next K subset in the ordering.
//    If the input T was the last in the ordering, then the output T
//    will be the first.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int ierror;
  int j;
// 
//  Return the first element.
// 
  if ( rank == - 1 )
  {
    i4vec_indicator ( k, t );
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = ksubset_lex_check ( k, n, t );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "KSUBSET_REVDOOR_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  j = 0;

  for ( ; ; )
  {
    if ( 0 < j || ( k % 2 ) == 0 )
    {
      j = j + 1;

      if ( k < j )
      {
        t[k-1] = k;
        rank = 0;
        return;
      }

      if ( t[j-1] != j )
      {
        t[j-1] = t[j-1] - 1;

        if ( j != 1 )
        {
          t[j-2] = j - 1;
        }
        rank = rank + 1;
        return;
      }
    }
    j = j + 1;

    if ( j < k )
    {
      if ( t[j-1] != t[j] - 1 )
      {
        break;
      }
    }
    else
    {
      if ( t[j-1] != n )
      {
        break;
      }
    }
  }

  t[j-1] = t[j-1] + 1;

  if ( j != 1 )
  {
    t[j-2] = t[j-1] - 1;
  }

  rank = rank + 1;

  return;
}
//****************************************************************************80

int *ksubset_revdoor_unrank ( int rank, int k, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    KSUBSET_REVDOOR_UNRANK computes the K subset of given revolving door rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the K subset.
// 
//    Input, int K, the number of elements each K subset must
//    have.  1 <= K <= N.
// 
//    Input, int N, the number of elements in the master set.
//    N must be positive.
// 
//    Output, int KSUBSET_REVDOOR_UNRANK[K], describes the K subset of the given
//    rank.  T(I) is the I-th element.  The elements must be listed in
//    ascending order.
// 
{
  int i;
  int nksub;
  int rank_copy;
  int *t;
  int x;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "KSUBSET_REVDOOR_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  if ( k < 1 || n < k )
  {
    cout << "\n";
    cout << "KSUBSET_REVDOOR_UNRANK - Fatal error!\n";
    cout << "  Input K is illegal.\n";
    exit ( 1 );
  }

  nksub = ksubset_enum ( k, n );

  if ( rank < 0 || nksub < rank )
  {
    cout << "\n";
    cout << "KSUBSET_REVDOOR_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  rank_copy = rank;

  t = new int[k];

  x = n;

  for ( i = k; 1 <= i; i-- )
  {
    while ( rank_copy < i4_choose ( x, i ) )
    {
      x = x - 1;
    }

    t[i-1] = x + 1;
    rank_copy = i4_choose ( x + 1, i ) - rank_copy - 1;
  }
  return t;
}
//****************************************************************************80

void marriage ( int n, int prefer[], int rank[], int fiancee[], int next[] )

//****************************************************************************80
// 
//  Purpose:
//
//    MARRIAGE finds a stable set of marriages for given preferences.
// 
//  Discussion:
// 
//    Given a set of N men and N women who must be married in pairs,
//    and information defining the relative rankings that each person
//    assigns to the candidates of the opposite sex, this routine finds
//    a stable set of marriages for them.
// 
//    A stable set of marriages is a pairing of the men and women with
//    the stability property: if M1 marries W1 and M2 marries W2, then
//    it is never the case that M1 and W2 would both prefer to be married
//    to each other.
// 
//    An important application of stable marriage algorithms occurs in
//    the annual matching of medical residents to hospitals.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Reference:
// 
//    Robert Sedgewick,
//    Algorithms in C,
//    Addison-Wesley, 1990,
//    ISBN: 0-201-51425-7,
//    LC: QA76.73.C15S43.
// 
//  Parameters:
// 
//    Input, int N, the number of pairs of men and women.
// 
//    Input, int PREFER[N*N]; for man I, the value of
//    PREFER(I,J) represents his J-th preference for a wife.
// 
//    Input, int RANK[N*N]; for woman I, the value of RANK(I,J)
//    represents her ranking of man number J.  A value of 1 for RANK(I,J)
//    means woman I ranks man J most preferable, while a value of N
//    would mean she ranked him least preferable.
// 
//    Output, int FIANCEE[N]; for woman I, FIANCEE(I) is the
//    man to whom she is now engaged.
// 
//    Output, int NEXT[N]; for man I, NEXT(I) is his preference
//    ranking for the woman to whom he is now engaged.  A value of 1 represents
//    his first choice, a value of N his last.
// 
{
  int i;
  int m;
  int temp;
  int w;
// 
//  For man I, NEXT(I) is the woman I has most recently proposed to,
//  and hence NEXT(I)+1 is the next one to try.
// 
  for ( i = 0; i < n; i++ )
  {
    next[i] = 0;
  }
// 
//  For woman I, FIANCEE(I) is the man she has agree to marry,
//  or 0 if she has not agreed to any man yet.
// 
  for ( i = 0; i < n; i++ )
  {
    fiancee[i] = -1;
  }
// 
//  Start with an unengaged man, and end with an engaged woman.
// 
  for ( i = 1; i <= n; i++ )
  {
    m = i;

    for ( ; ; )
    {
      next[m-1] = next[m-1] + 1;

      w = prefer[m-1+(next[m-1]-1)*n];

      if ( fiancee[w-1] == -1 )
      {
        fiancee[w-1] = m;
        break;
      }

      if ( rank[w-1+(m-1)*n] < rank[w-1+(fiancee[w-1]-1)*n] )
      {
        temp         = fiancee[w-1];
        fiancee[w-1] = m;
        m            = temp;
      }
    }
  }
  return;
}
//****************************************************************************80

int mountain ( int n, int x, int y )

//****************************************************************************80
// 
//  Purpose:
//
//    MOUNTAIN enumerates the mountains.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, ...
//    N must be positive.
// 
//    Input, int X, Y, ...
//    0 <= X <= 2 * N,
// 
//    Output, int MOUNTAIN, the value of the "mountain function"
//    M ( N, X, Y ), which is the number of all mountain ranges from
//    (X,Y) to (2*N,0) which do not drop below sea level.
// 
{
  int a;
  int b;
  int c;
  int value;
// 
//  Check.
// 
  if ( n <= 0 )
  {
    cout << "\n";
    cout << "MOUNTAIN - Fatal error!\n";
    cout << "  N <= 0.\n";
    cout << "  N = " << n << "\n";
    exit ( 1 );
  }
  else if ( x < 0 )
  {
    cout << "\n";
    cout << "MOUNTAIN - Fatal error!\n";
    cout << "  X < 0.\n";
    cout << "  X = " << x << "\n";
    exit ( 1 );
  }
  else if ( 2 * n < x )
  {
    cout << "\n";
    cout << "MOUNTAIN - Fatal error!\n";
    cout << "  2 * N < X.\n";
    cout << "  X = " << x << "\n";
    cout << "  N = " << n << "\n";
    exit ( 1 );
  }
// 
//  Special cases.
// 
  if ( y < 0 )
  {
    value = 0;
  }
  else if ( 2 * n < x + y )
  {
    value = 0;
  }
  else if ( ( ( x + y ) % 2 ) == 1 )
  {
   value = 0;
  }
  else
  {
    a = 2 * n - x;
    b = n - ( x + y ) / 2;
    c = n - 1 - ( x + y ) / 2;
    value = i4_choose ( a, b ) - i4_choose ( a, c );
  }
  return value;
}
//****************************************************************************80

int npart_enum ( int n, int npart )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_ENUM enumerates the number of partitions of N with NPART parts.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the integer to be partitioned.
//    Normally N must be positive, but for this routine any
//    N is allowed.
// 
//    Input, int NPART, the number of parts of the partition.
//    Normally, 1 <= NPART <= N is required,
//    but for this routine any value of NPART is allowed.
// 
//    Output, int NPART_ENUM is the number of partitions of N
//    with NPART parts.
// 
{
  int *p;
  int value;

  if ( n <= 0 )
  {
    value = 0;
  }  
  else if ( npart <= 0 || n < npart )
  {
    value = 0;
  }
  else
  {
    p = npart_table ( n, npart );

    value = p[n+npart*(n+1)];

    delete [] p;
  }

  return value;
}
//****************************************************************************80

int *npart_rsf_lex_random ( int n, int npart, int *seed )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_RSF_LEX_RANDOM returns a random RSF NPART partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input/output, int *SEED, a seed for the random number
//    generator.
// 
//    Output, int NPART_RSF_LEX_RANDOM[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
{
  int *a;
  int npartitions;
  int rank;

  npartitions = npart_enum ( n, npart );

  rank = i4_uniform ( 1, npartitions, seed );

  a = npart_rsf_lex_unrank ( rank, n, npart );

  return a;
}
//****************************************************************************80

int npart_rsf_lex_rank ( int n, int npart, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_RSF_LEX_RANK computes the lex rank of an RSF NPART partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
//    Output, int NPART_RSF_LEX_RANK, the rank of the partition.
// 
{
  int *b;
  int i;
  int ierror;
  int ncopy;
  int npartcopy;
  int *p;
  int rank;
// 
//  Check.
// 
  ierror = part_rsf_check ( n, npart, a );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "NPART_RSF_LEX_RANK - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Get the table of partitions of N with NPART parts.
// 
  p = npart_table ( n, npart );
// 
//  Copy the partition "backwards".
//
  b = new int[npart];

  for ( i = 1; i <= npart; i++ )
  {
    b[i-1] = a[npart-i];
  }

  rank = 0;
  ncopy = n;
  npartcopy = npart;

  while ( 0 < ncopy && 0 < npartcopy )
  {
    if ( b[npartcopy-1] == 1 )
    {
      ncopy = ncopy - 1;
      npartcopy = npartcopy - 1;
    }
    else
    {
      for ( i = 0; i < npartcopy; i++ )
      {
        b[i] = b[i] - 1;
      }
      rank = rank + p[ncopy-1+(npartcopy-1)*(n+1)];
      ncopy = ncopy - npartcopy;
    }
  }
  delete [] b;
  delete [] p;

  return rank;
}
//****************************************************************************80

void npart_rsf_lex_successor ( int n, int npart, int a[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_RSF_LEX_SUCCESSOR computes the RSF lex successor for NPART partitions.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be at least 1.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input/output, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int d;
  int i;
  int ierror;
  int j;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    if ( npart < 1 )
    {
      cout << "\n";
      cout << "NPART_RSF_LEX_SUCCESSOR - Fatal error!\n";
      cout << "  NPART < 1.\n";
      exit ( 1 );
    }

    for ( i = 0; i < npart - 1; i++ )
    {
      a[i] = 1;
    }
    a[npart-1] = n - ( npart - 1 );

    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = part_rsf_check ( n, npart, a );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "NPART_RSF_LEX_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Find the first index I for which A(NPART+1-I) + 1 < A(NPART).
// 
  i = 2;

  for ( ; ; )
  {
    if ( npart < i )
    {
      break;
    }

    if ( a[npart-i] + 1 < a[npart-1] )
    {
      break;
    }
    i = i + 1;
  }
// 
//  If no such index, we''ve reached the end of the line.
// 
  if ( i == npart + 1 )
  {
    for ( i = 0; i < npart - 1; i++ )
    {
      a[i] = 1;
    }
    a[npart-1] = n - ( npart - 1 );

    rank = 0;
    return;
  }
// 
//  Otherwise, increment A(NPART+1-I), and adjust other entries.
// 
  else
  {
    a[npart-i] = a[npart-i] + 1;
    d = - 1;

    for ( j = i - 1; 2 <= j; j-- )
    {
      d = d + a[npart-j] - a[npart-i];
      a[npart-j] = a[npart-i];
    }
    a[npart-1] = a[npart-1] + d;
  }
  rank = rank + 1;

  return;
}
//****************************************************************************80

int *npart_rsf_lex_unrank ( int rank, int n, int npart )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_RSF_LEX_UNRANK unranks an RSF NPART partition in the lex ordering.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int RANK, the rank of the partition.
// 
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Output, int NPART_RSF_LEX_UNRANK[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
{
  int *a;
  int i;
  int ncopy;
  int npartcopy;
  int npartitions;
  int *p;
  int rank_copy;
// 
//  Check.
// 
  if ( n <= 0 )
  {
    cout << "\n";
    cout << "NPART_RSF_LEX_UNRANK - Fatal error!\n";
    cout << "  The input N is illegal.\n";
    exit ( 1 );
  }

  if ( npart < 1 || n < npart )
  {
    cout << "\n";
    cout << "NPART_RSF_LEX_UNRANK - Fatal error!\n";
    cout << "  The input NPART is illegal.\n";
    exit ( 1 );
  }

  npartitions = npart_enum ( n, npart );

  if ( rank < 0 || npartitions < rank )
  {
    cout << "\n";
    cout << "NPART_RSF_LEX_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }
// 
//  Get the table of partitions of N with NPART parts.
// 
  p = npart_table ( n, npart );

  a = new int[npart];

  for ( i = 0; i < npart; i++ )
  {
    a[i] = 0;
  }

  rank_copy = rank;
  ncopy = n;
  npartcopy = npart;

  while ( 0 < ncopy )
  {
    if ( rank_copy < p[ncopy-1+(npartcopy-1)*(n+1)] )
    {
      a[npart-npartcopy] = a[npart-npartcopy] + 1;
      ncopy = ncopy - 1;
      npartcopy = npartcopy - 1;
    }
    else
    {
      for ( i = 1; i <= npartcopy; i++ )
      {
        a[npart-i] = a[npart-i] + 1;
      }
      rank_copy = rank_copy - p[ncopy-1+(npartcopy-1)*(n+1)];
      ncopy = ncopy - npartcopy;
    }
  }
  return a;
}
//****************************************************************************80

void npart_sf_lex_successor ( int n, int npart, int a[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_SF_LEX_SUCCESSOR computes SF NPART partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input/output, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.  The values in A must be in DESCENDING order.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int indx;
  int temp;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    i4vec_part2 ( n, npart, a );
    rank = 0;
    return;
  }
// 
//  Check.
// 
  part_sf_check ( n, npart, a );
// 
//  Find the last entry that is 2 or more.
// 
  for ( i = npart; 1 <= i; i-- )
  {
    if ( 1 < a[i-1] )
    {
      indx = i;
      break;
    }
  }
// 
//  As long as the last nonunit occurs after the first position,
//  have it donate 1 to the left.
// 
  if ( 1 < indx )
  {
    a[indx-1] = a[indx-1] - 1;
    a[indx-2] = a[indx-2] + 1;
    indx = indx - 1;

    for ( ; ; )
    {
      if ( indx <= 1 )
      {
        break;
      }

      if ( a[indx-1] <= a[indx-2] )
      {
        break;
      }

      temp      = a[indx-1];
      a[indx-1] = a[indx-2];
      a[indx-2] = temp;

      indx = indx - 1;
    }
// 
//  Sum the tail.
// 
    temp = 0;
    for ( i = indx; i < npart; i++ )
    {
      temp = temp + a[i];
    }
// 
//  Partition the tail sum equally over the tail.
// 
    i4vec_part2 ( temp, npart - indx, a+indx );

    rank = rank + 1;
  }
// 
//  If A(2) through A(NPART) are 1, then this is the last element.
//  Return the first one.
// 
  else
  {
    i4vec_part2 ( n, npart, a );
    rank = 0;
  }

  return;
}
//****************************************************************************80

int *npart_table ( int n, int npart )

//****************************************************************************80
// 
//  Purpose:
//
//    NPART_TABLE tabulates the number of partitions of N having NPART parts.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Output, int NPART_TABLE[(N+1)*(NPART+1)], P(I,J) is the number of
//    partitions of I having J parts.
// 
{
  int i;
  int j;
  int *p;

  p = new int[(n+1)*(npart+1)];

  p[0+0*(n+1)] = 1;
  for ( i = 1; i <= n; i++ )
  {
    p[i+0*(n+1)] = 0;
  }

  for ( i = 1; i <= n; i++ )
  {
    for ( j = 1; j <= npart; j++ )
    {
      if ( i < j )
      {
        p[i+j*(n+1)] = 0;
      }
      else if ( i < 2 * j )
      {
        p[i+j*(n+1)] = p[i-1+(j-1)*(n+1)];
      }
      else
      {
        p[i+j*(n+1)] = p[i-1+(j-1)*(n+1)] + p[i-j+j*(n+1)];
      }
    }
  }
  return p;
}
//****************************************************************************80

int part_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_ENUM enumerates the number of partitions of N.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the integer to be partitioned.
//    Normally N must be positive, but for this routine any
//    N is allowed.
// 
//    Output, int PART_ENUM is the number of partitions of N.
// 
{
  int *p;
  int value;

  if ( n < 0 )
  {
    value = 0;
  }
  else
  {
    p = part_table ( n );

    value = p[n];

    delete [] p;
  }
  return value;
}
//****************************************************************************80

int part_rsf_check ( int n, int npart, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_RSF_CHECK checks a reverse standard form partition of an integer.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.  The entries must be in ASCENDING order.
// 
//    Output, int PART_RSF_CHECK, error flag.
//    0, no error.
//    -1, N is illegal.
//    -2, NPART is illegal.
//    -3, the entries do not add up to N.
//    I, the I-th entry of A is illegal.
// 
{
  int i;
  int ierror;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  if ( npart < 1 || n < npart )
  {
    ierror = -2;
    return ierror;
  }
// 
//  Every entry must lie between 1 and N.
// 
  for ( i = 0; i < npart; i++ )
  {
    if ( a[i] < 1 || n < a[i] )
    {
      ierror = i + 1;
      return ierror;
    }
  }
// 
//  The entries must be in ascending order.
// 
  for ( i = 1; i < npart; i++ )
  {
    if ( a[i] < a[i-1] )
    {
      ierror = i + 1;
      return ierror;
    }
  }
// 
//  The entries must add up to N.
// 
  if ( i4vec_sum ( npart, a ) != n )
  {
    ierror = -3;
    return ierror;
  }

  return ierror;
}
//****************************************************************************80

void part_sf_check ( int n, int npart, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_SF_CHECK checks a standard form partition of an integer.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    27 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.  The entries must be in DESCENDING order.
//
{
  int i;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "PART_SF_CHECK - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  if ( npart < 1 || n < npart )
  {
    cerr << "\n";
    cerr << "PART_SF_CHECK - Fatal error!\n";
    cerr << "  1 <= NPART <= N is required,\n";
    cerr << "  but NPART = " << npart << "\n";
    exit ( 1 );
  }
// 
//  Every entry must lie between 1 and N.
// 
  for ( i = 0; i < npart; i++ )
  {
    if ( a[i] < 1 || n < a[i] )
    {
      cerr << "\n";
      cerr << "PART_SF_CHECK - Fatal error!\n";
      cerr << "  1 <= A[*] <= N is required.\n";
      cerr << "  but A[" << i << "] = " << a[i] << "\n";
      exit ( 1 );
    }
  }
// 
//  The entries must be in descending order.
// 
  for ( i = 1; i < npart; i++ )
  {
    if ( a[i-1] < a[i] )
    {
      cerr << "\n";
      cerr << "PART_SF_CHECK - Fatal error!\n";
      cerr << "  The entries are not in descending order.\n";
      exit ( 1 );
    }
  }
// 
//  The entries must add up to N.
// 
  if ( i4vec_sum ( npart, a ) != n )
  {
    cerr << "\n";
    cerr << "PART_SF_CHECK - Fatal error!\n";
    cerr << "  The NPART entries of A don't add up to N.\n";
    cerr << "  N = " << n << "\n";
    cerr << "  NPART = " << npart << "\n";
    i4vec_transpose_print ( npart, a, "  A:" );
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int *part_sf_conjugate ( int n, int npart, int a[], int &npart2 )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_SF_CONJUGATE computes the conjugate of a partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input, int A[N], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
//    Output, int &NPART2, the number of parts of the conjugate
//    partition.
// 
//    Output, int PART_SF_CONJUGATE[N], contains the conjugate partition.
// 
{
  int *b;
  int i;
  int j;
// 
//  Check.
// 
  part_sf_check ( n, npart, a );

  npart2 = a[0];

  b = new int[n];

  for ( i = 0; i < npart2; i++ )
  {
    b[i] = 0;
  }

  for ( i = 0; i < npart; i++ )
  {
    for ( j = 0; j < a[i]; j++ )
    {
      b[j] = b[j] + 1;
    }
  }

  return b;
}
//****************************************************************************80

int part_sf_majorize ( int n, int nparta, int a[], int npartb, int b[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_SF_MAJORIZE determines if partition A majorizes partition B.
// 
//  Discussion:
// 
//    The partitions must be in standard form.
// 
//    If A, with NPARTA parts, and B, with NPARTB parts, are both partitions
//    of the same positive integer N, then we say that A majorizes B if,
//    for every index K from 1 to N, it is true that
// 
//      sum ( 1 <= I <= K ) B(I) <= sum ( 1 <= I <= K ) A(I)
// 
//    where entries of A beyond index NPARTA, and of B beyond BPARTB
//    are assumed to be 0.  We say that A strictly majorizes B if
//    A majorizes B, and for at least one index K the inequality is strict.
// 
//    For any two partitions of N, it is possible that A majorizes B,
//    B majorizes A, both partitions majorize each other (in which case
//    they are equal), or that neither majorizes the other.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Reference:
// 
//    Jack vanLint, Richard Wilson,
//    A Course in Combinatorics,
//    Cambridge, 1992,
//    ISBN: 0-521-42260-4,
//    LC: QA164.L56.
// 
//  Parameters:
// 
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NPARTA, the number of parts in partition A.
//    1 <= NPARTA <= N.
// 
//    Input, int A[NPARTA], contains partition A in standard
//    form.  A(1) through A(NPARTA) contain nonzero integers which sum to N.
// 
//    Input, int NPARTB, the number of parts in partition B.
//    1 <= NPARTB <= N.
// 
//    Input, int B[NPARTB], contains partition B in standard
//    form.  B(1) through B(NPARTB) contain nonzero integers which sum to N.
// 
//    Output, int PART_SF_MAJORIZE, the result of the comparison.
//    -2, A and B are incomparable, but would have been -1.
//    -1, A < B, (A is strictly majorized by B),
//     0, A = B, (A and B are identical),
//    +1, A > B, (A strictly majorizes B),
//    +2, A and B are incomparable, but would have been +1.
// 
{
  int i;
  int result;
  int suma;
  int sumb;
// 
//  Check.
// 
  part_sf_check ( n, nparta, a );

  part_sf_check ( n, npartb, b );

  result = 0;
  suma = 0;
  sumb = 0;

  for ( i = 0; i < i4_min ( nparta, npartb ); i++ )
  {
    if ( i < nparta )
    {
      suma = suma + a[i];
    }

    if ( i < npartb )
    {
      sumb = sumb + b[i];
    }

    if ( result == -1 )
    {
      if ( sumb < suma )
      {
        result = -2;
        return result;
      }
    }
    else if ( result == 0 )
    {
      if ( suma < sumb )
      {
        result = -1;
      }
      else if ( sumb < suma )
      {
        result = +1;
      }
    }
    else if ( result == + 1 )
    {
      if ( suma < sumb )
      {
        result = +2;
        return result;
      }
    }
  }

  return result;
}
//****************************************************************************80

void part_successor ( int n, int &npart, int a[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_SUCCESSOR computes the lexicographic partition successor.
// 
//  Discussion:
// 
//    PART_SUCCESSOR is "inspired by" the GenPartitions algorithm,
//    but instead of relying on recursion, generates the partitions
//    one at a time.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input/output, int &NPART, the number of parts of the
//    partition.  1 <= NPART <= N.
// 
//    Input/output, int A[N], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int asum;
  int i;
  int ihi;
  int j;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 1;
    }
    npart = n;
    rank = 0;
    return;
  }
// 
//  Check.
// 
  part_sf_check ( n, npart, a );
// 
//  If possible, increment the first intermediate position that
//  is less than its left hand neighbor, and has at least one
//  right hand neighbor.
// 
  ihi = npart - 1;

  for ( i = ihi; 2 <= i; i-- )
  {
    if ( a[i-1] < a[i-2] )
    {
      asum = - 1;
      for ( j = i + 1; j <= npart; j++ )
      {
        asum = asum + a[j-1];
      }
      a[i-1] = a[i-1] + 1;
      for ( j = i + 1; j <= npart; j++ )
      {
        a[j-1] = 0;
      }
      npart = i + asum;
      for ( j = i + 1; j <= npart; j++ )
      {
        a[j-1] = 1;
      }
      rank = rank + 1;
      return;
    }
  }
// 
//  A) there are two or more parts
//  Increment the first, replace the rest by 1''s.
// 
  if ( 2 <= npart )
  {
    a[0] = a[0] + 1;
    for ( j = 2; j <= npart; j++ )
    {
      a[j-1] = 0;
    }
    npart = n - a[0] + 1;
    for ( j = 2; j <= npart; j++ )
    {
      a[j-1] = 1;
    }
    rank = rank + 1;
  }
// 
//  B) there is only one part.
//  We have reached the last item.
//  Return the first one.
// 
  else if ( npart == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 1;
    }
    npart = n;
    rank = 0;
  }
  return;
}
//****************************************************************************80

int *part_table ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PART_TABLE tabulates the number of partitions of N.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Output, int P[N+1], P(I) is the number of partitions of I.
// 
{
  int i;
  int j;
  int *p;
  int psum;
  int sign;
  int w;
  int wprime;

  p = new int[n+1];

  p[0] = 1;
  p[1] = 1;

  for ( i = 2; i <= n; i++ )
  {
    sign = 1;
    psum = 0;
    w = 1;
    j = 1;
    wprime = w + j;

    while ( w < n )
    {
      if ( 0 <= i - w )
      {
        if ( sign == 1 )
        {
          psum = psum + p[i-w];
        }
        else
        {
          psum = psum - p[i-w];
        }
      }

      if ( wprime <= i )
      {
        if ( sign == 1 )
        {
          psum = psum + p[i-wprime];
        }
        else
        {
          psum = psum - p[i-wprime];
        }
      }
      w = w + 3 * j + 1;
      j = j + 1;
      wprime = w + j;
      sign = - sign;
    }
    p[i] = psum;
  }
  return p;
}
//****************************************************************************80

int *partition_greedy ( int n, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PARTITION_GREEDY attacks the partition problem with a greedy algorithm.
// 
//  Discussion:
// 
//    Given a collection of N not necessarily distinct positive integers A(I),
//    it is desired to partition the values into two groups, whose sums are
//    as close as possible.
// 
//  Algorithm:
// 
//    Begin with sets 1 and 2 empty.
// 
//    Process the data in descending order of magnitude.
// 
//    The next item A(I) is added to set 1 or set 2, whichever has the
//    smallest current sum.
// 
//    Stop as soon as all items have been allocated.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Reference:
// 
//    Brian Hayes,
//    The Easiest Hard Problem,
//    American Scientist,
//    Volume 90, Number 2, March-April 2002, pages 113-117.
// 
//  Parameters:
// 
//    Input, int N, the number of values.  N must be positive.
// 
//    Input/output, int A[N], a collection of positive values.
//    On output, A has been sorted into descending order.
// 
//    Output, int PARTITION_GREEDY[N]; an entry is 0 if A[I] is part of
//    set 0, and 1 if it is assigned to set 1.
// 
{
  int i;
  int *indx;
  int j;
  int sums[2];

  sums[0] = 0;
  sums[1] = 0;

  i4vec_sort_insert_d ( n, a );

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    if ( sums[0] < sums[1] )
    {
      j = 0;
    }
    else
    {
      j = 1;
    }
    indx[i] = j;
    sums[j] = sums[j] + a[i];
  }
  return indx;
}
//****************************************************************************80

int partn_enum ( int n, int nmax )

//****************************************************************************80
// 
//  Purpose:
//
//    PARTN_ENUM enumerates the partitions of N with maximum element NMAX.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the integer to be partitioned.
//    Normally N must be positive, but for this routine any
//    N is allowed.
// 
//    Input, int NMAX, the maximum element in the partition.
//    Normally, 1 <= NMAX <= N is required,
//    but for this routine any value of NMAX is allowed.
// 
//    Output, int PARTN_ENUM is the number of partitions of N
//    with maximum element NMAX.
// 
{
  int *p;
  int value;

  if ( n <= 0 )
  {
    value = 0;
  }
  else if ( nmax <= 0 || n < nmax )
  {
    value = 0;
  }
  else
  {
    p = npart_table ( n, nmax );

    value = p[n+nmax*(n+1)];

    delete [] p;
  }

  return value;
}
//****************************************************************************80

int partn_sf_check ( int n, int nmax, int npart, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PARTN_SF_CHECK checks an SF partition of an integer with largest entry NMAX.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NMAX, the value of the largest entry.
//    1 <= NMAX <= N.
// 
//    Input, int NPART, the number of parts of the partition.
//    1 <= NPART <= N.
// 
//    Input, int A[NPART], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.  The entries must be in DESCENDING order.
// 
//    Output, int PARTN_SF_CHECK, error flag.
//    0, no error.
//    -1, N is illegal.
//    -2, NMAX is illegal.
//    -3, NPART is illegal.
//    -3, the entries do not add up to N.
//    I, the I-th entry of A is illegal.
// 
{
  int asum;
  int i;
  int ierror;

  ierror = 0;

  if ( n < 1 )
  {
    ierror = -1;
    return ierror;
  }

  if ( nmax < 1 || n < nmax )
  {
    ierror = -2;
    return ierror;
  }

  if ( npart < 1 || n < npart )
  {
    ierror = -3;
    return ierror;
  }
// 
//  Entry 1 must be NMAX.
// 
  if ( a[0] != nmax )
  {
    ierror = 1;
    return ierror;
  }
// 
//  Every entry must lie between 1 and N.
// 
  for ( i = 0; i < npart; i++ )
  {
    if ( a[i] < 1 || n < a[i] )
    {
      ierror = i + 1;
      return ierror;
    }
  }
// 
//  The entries must be in descending order.
// 
  for ( i = 1; i < npart; i++ )
  {
    if ( a[i-1] < a[i] )
    {
      ierror = i + 1;
      return ierror;
    }
  }
// 
//  The entries must add up to N.
// 
  asum = i4vec_sum ( npart, a );

  if ( asum != n )
  {
    ierror = -3;
    return ierror;
  }

  return ierror;
}
//****************************************************************************80

void partn_successor ( int n, int nmax, int &npart, int a[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    PARTN_SUCCESSOR computes partitions whose largest part is NMAX.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int N, the integer to be partitioned.
//    N must be positive.
// 
//    Input, int NMAX, the maximum size of any part of the
//    partition.  1 <= NMAX <= N.
// 
//    Input/output, int NPART, the number of parts of the
//    partition.  1 <= NPART <= N.
// 
//    Input/output, int A[N], contains the partition.
//    A(1) through A(NPART) contain the nonzero integers which
//    sum to N.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
  int index;
  int temp;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    npart = n + 1 - nmax;
    a[0] = nmax;
    for ( i = 1; i < npart; i++ )
    {
      a[i] = 1;
    }
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = partn_sf_check ( n, nmax, npart, a );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "PARTN_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  If there are at least two parts, and the next to last is not NMAX,
//  then rob the last part and pay the next to the last part.
//  Then, if the next to last part is too big, swap it leftwards.
// 
  if ( 1 < npart )
  {
    if ( a[npart-2] < nmax )
    {
      a[npart-1] = a[npart-1] - 1;
      a[npart-2] = a[npart-2] + 1;
      index = npart - 1;

      for ( ; ; )
      {
        if ( index <= 1 )
        {
          break;
        }

        if ( a[index-1] <= a[index-2] )
        {
          break;
        }

        temp       = a[index-2];
        a[index-2] = a[index-1];
        a[index-1] = temp;

        index = index - 1;
      }
// 
//  Sum the tail.
// 
      temp = 0;
      for ( i = index; i < npart; i++ )
      {
        temp = temp + a[i];
      }
// 
//  Spread the sum as 1''s.
// 
      npart = index + temp;
      for ( i = index; i < npart; i++ )
      {
        a[i] = 1;
      }
      rank = rank + 1;
      return;
    }
  }
// 
//  Otherwise, we have reached the last item.
//  Return the first one.
// 
  else
  {
    npart = n + 1 - nmax;
    a[0] = nmax;
    for ( i = 1; i < npart; i++ )
    {
      a[i] = 1;
    }
    rank = 0;
    return;
  }

  return;
}
//****************************************************************************80

void perm_check ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_CHECK checks a representation of a permutation.
// 
//  Discussion:
// 
//    The routine is given N and P, a vector of length N.
//    P is a legal represention of a permutation of the integers from
//    1 to N if and only if every integer from 1 to N occurs
//    as a value of P(I) for some I between 1 and N.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], the array to check.
// 
{
  int i;
  int ifind;
  int iseek;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "PERM_CHECK - Fatal error!\n";
    cerr << "  Input N = " << n << " < 1.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( p[i] < 1 || n < p[i] )
    {
      cerr << "\n";
      cerr << "PERM_CHECK - Fatal error!\n";
      cerr << "  P[" << i << "] = " << p[i] << "\n";
      cerr << "  but 1 <= p[i] <= " << n << " is required.\n";
      exit ( 1 );
    }
  }

  for ( iseek = 1; iseek <= n; iseek++ )
  {
    ifind = -1;
    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == iseek )
      {
        ifind = i;
        break;
      }
    }

    if ( ifind == -1 )
    {
      cerr << "\n";
      cerr << "PERM_CHECK - Fatal error!\n";
      cerr << "  Could not locate the value " << iseek << "\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************80

int perm_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_ENUM enumerates the permutations on N digits.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of values being permuted.
//    N must be nonnegative.
// 
//    Output, int PERM_ENUM, the number of distinct elements.
// 
{
  int value;

  value = i4_factorial ( n );

  return value;
}
//****************************************************************************80

int *perm_inv ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_INV computes the inverse of a permutation.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], describes the permutation.
//    P(I) is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Output, int PERM_INV[N], the inverse permutation.
// 
{
  int i;
  int *pinv;
// 
//  Check.
// 
  perm_check ( n, p );

  pinv = new int[n];

  for ( i = 0; i < n; i++ )
  {
    pinv[p[i]-1] = i + 1;
  }

  return pinv;
}
//****************************************************************************80

int perm_lex_rank ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_LEX_RANK computes the lexicographic rank of a permutation.
// 
//  Discussion:
// 
//    The original code altered the input permutation.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], describes the permutation.
//    P[I] is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Output, int PERM_LEX_RANK, the rank of the permutation.
// 
{
  int i;
  int j;
  int *pcopy;
  int rank;
// 
//  Check.
// 
  perm_check ( n, p );

  rank = 0;
  pcopy = new int[n];

  for ( i = 0; i < n; i++ )
  {
    pcopy[i] = p[i];
  }

  for ( j = 0; j < n; j++ )
  {
    rank = rank + ( pcopy[j] - 1 ) * i4_factorial ( n - 1 - j );
    for ( i = j + 1; i < n; i++ )
    {
      if ( pcopy[j] < pcopy[i] )
      {
        pcopy[i] = pcopy[i] - 1;
      }
    }
  }
  delete [] pcopy;

  return rank;
}
//****************************************************************************80

void perm_lex_successor ( int n, int p[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_LEX_SUCCESSOR computes the lexicographic permutation successor.
// 
//  Example:
// 
//    RANK  Permutation
// 
//       0  1 2 3 4
//       1  1 2 4 3
//       2  1 3 2 4
//       3  1 3 4 2
//       4  1 4 2 3
//       5  1 4 3 2
//       6  2 1 3 4
//       ...
//      23  4 3 2 1
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input/output, int P[N], describes the permutation.
//    P(I) is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int j;
  int temp;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    i4vec_indicator ( n, p );
    rank = 0;
    return;
  }
// 
//  Check.
// 
  perm_check ( n, p );
// 
//  Seek I, the highest index for which the next element is bigger.
// 
  i = n - 1;

  for ( ; ; )
  {
    if ( i <= 0 )
    {
      break;
    }

    if ( p[i-1] <= p[i] )
    {
      break;
    }
    i = i - 1;
  }
// 
//  If no I could be found, then we have reach the final permutation,
//  N, N-1, ..., 2, 1.  Time to start over again.
// 
  if ( i == 0 )
  {
    i4vec_indicator ( n, p );
    rank = 0;
  }
  else
  {
// 
//  Otherwise, look for the the highest index after I whose element
//  is bigger than I''s.  We know that I+1 is one such value, so the
//  loop will never fail.
// 
    j = n;
    while ( p[j-1] < p[i-1] )
    {
      j = j - 1;
    }
// 
//  Interchange elements I and J.
// 
    temp = p[i-1];
    p[i-1] = p[j-1];
    p[j-1] = temp;
// 
//  Reverse the elements from I+1 to N.
// 
    i4vec_reverse ( n - i, p+i );

    rank = rank + 1;
  }

  return;
}
//****************************************************************************80

int *perm_lex_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_LEX_UNRANK computes the permutation of given lexicographic rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int RANK, the rank of the permutation.
// 
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Output, int PERM_LEX_UNRANK[N], describes the permutation.
// 
{
  int d;
  int i;
  int j;
  int nperm;
  int *p;
  int rank_copy;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "PERM_LEX_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  nperm = perm_enum ( n );

  if ( rank < 0 || nperm < rank )
  {
    cout << "\n";
    cout << "PERM_LEX_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  rank_copy = rank;

  p = new int[n];

  p[n-1] = 1;

  for ( j = 1; j <= n - 1; j++ )
  {
    d = ( rank_copy % i4_factorial ( j + 1 ) ) / i4_factorial ( j );
    rank_copy = rank_copy - d * i4_factorial ( j );
    p[n-j-1] = d + 1;

    for ( i = n - j + 1; i <= n; i++ )
    {
      if ( d < p[i-1] )
      {
        p[i-1] = p[i-1] + 1;
      }
    }
  }
  return p;
}
//****************************************************************************80

int *perm_mul ( int n, int p[], int q[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_MUL computes the product of two permutations.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Reference:
// 
//    Donald Kreher, Douglas Simpson,inson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
// 
//  Parameters:
// 
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], Q[N], describes the permutation factors.
// 
//    Output, int PERMN_MUL[N], the product permutation P * Q.
//    R(I) = P(Q(I)).
// 
{
  int i;
  int *r;
// 
//  Check.
// 
  perm_check ( n, p );

  perm_check ( n, q );
// 
//  Use a temporary vector for the result, to avoid problems if
//  some arguments are actually identified.
// 
  r = new int[n];

  for ( i = 0; i < n; i++ )
  {
    r[i] = p[q[i]-1];
  }

  return r;
}
//****************************************************************************80

int perm_parity ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_PARITY computes the parity of a permutation.
// 
//  Discussion:
// 
//    The routine requires the use of a temporary array.
// 
//    A permutation is called "even" or "odd", depending on whether
//    it is equivalent to an even or odd number of pairwise
//    transpositions.  This is known as the "parity" of the
//    permutation.
// 
//    The "sign" of a permutation is +1 if it has even parity,
//    and -1 if it has odd parity.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], describes the permutation.
//    P(I) is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Output, int PERM_PARITY, the parity of the permutation.
//    0, the permutation has even parity.
//    1, the permutation has odd parity.
// 
{
  int *a;
  int c;
  int i;
  int j;
  int parity;
// 
//  Check.
// 
  perm_check ( n, p );

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }

  c = 0;

  for ( j = 1; j <= n; j++ )
  {
    if ( a[j-1] == 0 )
    {
      c = c + 1;
      a[j-1] = 1;
      i = j;

      while ( p[i-1] != j )
      {
        i = p[i-1];
        a[i-1] = 1;
      }
    }
  }

  parity = ( n - c ) % 2;

  delete [] a;

  return parity;
}

//****************************************************************************80

void perm_print ( int n, int p[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_PRINT prints a permutation.
//
//  Discussion:
//
//    The permutation is assumed to be zero-based.
//
//  Example:
//
//    Input:
//
//      P = 6 1 2 0 4 2 5
//
//    Printed output:
//
//      "This is the permutation:"
//
//      0 1 2 3 4 5 6
//      6 1 2 0 4 2 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects permuted.
//
//    Input, int P[N], the permutation, in standard index form.
//
//    Input, string TITLE, a title.
//    If no title is supplied, then only the permutation is printed.
//
{
  int i;
  int ihi;
  int ilo;
  int inc = 20;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";

    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      cout << "\n";
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << i;
      }
      cout << "\n";
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << p[i];
      }
      cout << "\n";
    }
  }
  else
  {
    for ( ilo = 0; ilo < n; ilo = ilo + inc )
    {
      ihi = ilo + inc;
      if ( n < ihi ) 
      {
        ihi = n;
      }
      cout << "  ";
      for ( i = ilo; i < ihi; i++ )
      {
        cout << setw(4) << p[i];
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

int perm_tj_rank ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_TJ_RANK computes the Trotter-Johnson rank of a permutation.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], describes the permutation.
//    P(I) is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Output, int PERM_TJ_RANK, the rank of the permutation.
// 
{
  int i;
  int j;
  int k;
  int rank;
// 
//  Check.
// 
  perm_check ( n, p );

  rank = 0;

  for ( j = 2; j <= n; j++ )
  {
    k = 1;
    i = 1;

    while ( p[i-1] != j )
    {
      if ( p[i-1] < j )
      {
        k = k + 1;
      }
      i = i + 1;
    }

    if ( ( rank % 2 ) == 0 )
    {
      rank = j * rank + j - k;
    }
    else
    {
      rank = j * rank + k - 1;
    }
  }

  return rank;
}
//****************************************************************************80

void perm_tj_successor ( int n, int p[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_TJ_SUCCESSOR computes the Trotter-Johnson permutation successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input/output, int P[N], describes the permutation.
//    P(I) is the item which is permuted into the I-th place
//    by the permutation.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int d;
  bool done;
  int i;
  int m;
  int par;
  int *q;
  int st;
  int temp;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    i4vec_indicator ( n, p );
    rank = 0;
    return;
  }
// 
//  Check.
// 
  perm_check ( n, p );

  q = new int[n];

  st = 0;
  for ( i = 0; i < n; i++ )
  {
    q[i] = p[i];
  }
  done = false;
  m = n;

  while ( 1 < m && !done )
  {
    d = 1;
    while ( q[d-1] != m )
    {
      d = d + 1;
    }

    for ( i = d; i < m; i++ )
    {
      q[i-1] = q[i];
    }

    par = perm_parity ( m - 1, q );

    if ( par == 1 )
    {
      if ( d == m )
      {
        m = m - 1;
      }
      else
      {
        temp      = p[st+d-1];
        p[st+d-1] = p[st+d];
        p[st+d]   = temp;
        done = true;
      }
    }
    else
    {
      if ( d == 1 )
      {
        m = m - 1;
        st = st + 1;
      }
      else
      {
        temp      = p[st+d-1];
        p[st+d-1] = p[st+d-2];
        p[st+d-2] = temp;
        done = true;
      }
    }
  }
// 
//  Last element was input.  Return first one.
// 
  if ( m == 1 )
  {
    i4vec_indicator ( n, p );
    rank = 0;
    delete [] q;
    return;
  }

  rank = rank + 1;

  delete [] q;

  return;
}
//****************************************************************************80

int *perm_tj_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_TJ_UNRANK computes the permutation of given Trotter-Johnson rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the permutation.
// 
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Output, int PERM_TJ_UNRANK[N], describes the permutation.
// 
{
  int i;
  int j;
  int k;
  int jhi;
  int nperm;
  int *p;
  int r1;
  int r2;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "PERM_TJ_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  nperm = perm_enum ( n );

  if ( rank < 0 || nperm < rank )
  {
    cout << "\n";
    cout << "PERM_TJ_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  p = new int[n];

  p[0] = 1;
  r2 = 0;

  for ( j = 2; j <= n; j++ )
  {
// 
//  Replace this ratio of factorials!
// 
    r1 = ( rank * i4_factorial ( j ) ) / i4_factorial ( n );
    k = r1 - j * r2;

    if ( ( r2 % 2 ) == 0 )
    {
      jhi = j - k;
    }
    else
    {
      jhi = k + 1;
    }

    for ( i = j - 1; jhi <= i; i-- )
    {
      p[i] = p[i-1];
    }
    p[jhi-1] = j;

    r2 = r1;
  }

  return p;
}
//****************************************************************************80

void perm_to_cycle ( int n, int p[], int &ncycle, int t[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PERM_TO_CYCLE converts a permutation from array to cycle form.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int N, the number of values being permuted.
//    N must be positive.
// 
//    Input, int P[N], describes the permutation using a
//    single array.  For each index I, I -> P(I).
// 
//    Output, int &NCYCLE, the number of cycles.
//    1 <= NCYCLE <= N.
// 
//    Output, int T[N], INDEX[N], describes the permutation
//    as a collection of NCYCLE cycles.  The first cycle is
//    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
// 
{
  int i;
  int j;
  int nset;
// 
//  Check.
// 
  perm_check ( n, p );
// 
//  Initialize.
// 
  ncycle = 0;
  for ( i = 0; i < n; i++ )
  {
    index[i] = 0;
  }
  for ( i = 0; i < n; i++ )
  {
    t[i] = 0;
  }
  nset = 0;
// 
//  Find the next unused entry.
//
  for ( i = 1; i <= n; i++ )
  {
    if ( 0 < p[i-1] )
    {
      ncycle = ncycle + 1;
      index[ncycle-1] = 1;

      nset = nset + 1;
      t[nset-1] = p[i-1];
      p[i-1] = - p[i-1];

      for ( ; ; )
      {
        j = t[nset-1];

        if ( p[j-1] < 0 )
        {
          break;
        }

        index[ncycle-1] = index[ncycle-1] + 1;

        nset = nset + 1;
        t[nset-1] = p[j-1];
        p[j-1] = - p[j-1];
      }
    }
  }
// 
//  If no unused entries remain, we are done.
//  Restore the sign of the permutation and return.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }

  return;
}
//****************************************************************************80

int pruefer_check ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_CHECK checks a Pruefer code.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input, int P[N-2], the Pruefer code for the tree.
// 
//    Output, int PRUEFER_CHECK, error flag.
//    0, no error.
//    -1, N is less than 3.
//    J, the element P(J) is illegal.
// 
{
  int i;
  int ierror;

  ierror = 0;

  if ( n < 3 )
  {
    ierror = -1;
    return ierror;
  }

  for ( i = 0; i < n - 2; i++ )
  {
    if ( p[i] < 1 || n < p[i] )
    {
      ierror = i + 1;
      return ierror;
    }
  }
  return ierror;
}
//****************************************************************************80

int pruefer_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_ENUM enumerates the Pruefer codes on N-2 digits.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of digits in the code, plus 2.
//    N must be at least 3.
// 
//    Output, int NCODE, the number of distinct elements.
// 
{
  int value;

  value = i4_power ( n, n - 2 );

  return value;
}
//****************************************************************************80

int pruefer_rank ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_RANK ranks a Pruefer code.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input, int P[N-2], the Pruefer code for the tree.
// 
//    Output, int PRUEFER_RANK, the rank of the Pruefer code.
// 
{
  int i;
  int ierror;
  int k;
  int rank;
// 
//  Check.
// 
  ierror = pruefer_check ( n, p );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "PRUEFER_RANK - Fatal error!\n";
    cout << "  Input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  rank = 0;
  k = 1;
  for ( i = n - 3; 0 <= i; i-- )
  {
    rank = rank + k * ( p[i] - 1 );
    k = k * n;
  }

  return rank;
}
//****************************************************************************80

void pruefer_successor ( int n, int p[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_SUCCESSOR computes the lexical Pruefer sequence successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    27 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input/output, int P(N-2), on input, the Pruefer code
//    for a tree, and on output, its lexical successor.
// 
//    Input/output, int RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int i;
  int ierror;
  int j;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    for ( i = 0; i < n - 2; i++ )
    {
      p[i] = 1;
    }
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = pruefer_check ( n, p );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "PRUEFER_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }

  j = n - 2;

  for ( ; ; )
  {
    if ( p[j-1] != n )
    {
      break;
    }

    j = j - 1;

    if ( j <= 0 )
    {
      break;
    }
  }

  if ( j != 0 )
  {
    p[j-1] = p[j-1] + 1;
    for ( i = j + 1; i <= n - 2; i++ )
    {
      p[i-1] = 1;
    }
    rank = rank + 1;
  }
  else
  {
    for ( i = 0; i < n - 2; i++ )
    {
      p[i] = 1;
    }
    rank = 0;
  }
  return;
}
//****************************************************************************80

void pruefer_to_tree ( int n, int p[], int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_TO_TREE converts a Pruefer code to a tree.
// 
//  Discussion:
// 
//    The original code attempts to tack on an extra entry to P.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input, int P[N-2], the Pruefer code for the tree.
// 
//    Output, int T[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
{
  int *d;
  int i;
  int ierror;
  int j;
  int x;
  int y;
// 
//  Check.
// 
  ierror = pruefer_check ( n, p );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "PRUEFER_TO_TREE - Fatal error!\n";
    cout << "  The input array is illegal!\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Initialize the tree to 0.
// 
  for ( j = 0; j < n - 1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      t[i+j*2] = 0;
    }
  }

  d = new int[n];

  for ( i = 0; i < n; i++ )
  {
    d[i] = 1;
  }

  for ( i = 0; i < n - 2; i++ )
  {
    d[p[i]-1] = d[p[i]-1] + 1;
  }

  for ( i = 1; i <= n - 1; i++ )
  {
    x = n;
    while ( d[x-1] != 1 )
    {
      x = x - 1;
    }

    if ( i == n - 1 )
    {
      y = 1;
    }
    else
    {
      y = p[i-1];
    }

    d[x-1] = d[x-1] - 1;
    d[y-1] = d[y-1] - 1;

    t[0+(i-1)*2] = x;
    t[1+(i-1)*2] = y;
  }

  delete [] d;

  return;
}
//****************************************************************************80

int *pruefer_to_tree_new ( int n, int p[] )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_TO_TREE_NEW converts a Pruefer code to a tree.
// 
//  Discussion:
// 
//    The original code attempts to tack on an extra entry to P.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input, int P[N-2], the Pruefer code for the tree.
// 
//    Output, int PRUEFER_TO_TREE_NEW[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
{
  int *t;

  t = new int[2*(n-1)];

  pruefer_to_tree ( n, p, t );

  return t;
}
//****************************************************************************80

int *pruefer_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    PRUEFER_UNRANK unranks a Pruefer code.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int RANK, the rank of the Pruefer code.
// 
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Output, int P[N-2], the Pruefer code for the tree.
// 
{
  int i;
  int ncode;
  int *p;
  int rank_copy;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "PRUEFER_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  if ( n < 3 )
  {
    p = NULL;
    return p;
  }

  ncode = pruefer_enum ( n );

  if ( rank < 0 || ncode < rank )
  {
    cout << "\n";
    cout << "PRUEFER_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  rank_copy = rank;
  p = new int[n-2];

  for ( i = n - 3; 0 <= i; i-- )
  {
    p[i] = ( rank_copy % n ) + 1;
    rank_copy = ( rank_copy - p[i] + 1 ) / n;
  }
  return p;
}
//****************************************************************************80

void queens ( int n, int iarray[], int k, int &nstack, int istack[], 
  int maxstack )

//****************************************************************************80
// 
//  Purpose:
//
//    QUEENS finds possible positions for the K-th nonattacking queen.
// 
//  Discussion:
// 
//    The chessboard is N by N, and is being filled one column at a time,
//    with a tentative solution to the nonattacking queen problem.  So
//    far, K-1 rows have been chosen, and we now need to provide a list
//    of all possible rows that might be used in column K.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the total number of queens to place, and
//    the length of a side of the chessboard.
// 
//    Input, int IARRAY[N].  The first K-1 entries of IARRAY
//    record the rows into which queens have already been placed.
// 
//    Input, int K, the column for which we need possible
//    row positions for the next queen.
// 
//    Input/output, int &NSTACK, the current length of stack.
//    On output, this has been updated.
// 
//    Input/output, int ISTACK[MAXSTACK].  On output, we
//    have added the candidates, and the number of candidates, to the end
//    of the stack.
// 
//    Input, int MAXSTACK, maximum dimension of ISTACK.
// 
{
  bool diag;
  int irow;
  int jcol;
  int ncan;
  bool row;

  ncan = 0;

  for ( irow = 1; irow <= n; irow++ )
  {
// 
//  If row IROW has already been used, that is it.
// 
    row = false;

    for ( jcol = 1; jcol <= k - 1; jcol++ )
    {
      if ( iarray[jcol-1] == irow )
      {
        row = true;
      }
    }

    if ( !row )
    {
      diag = false;

      for ( jcol = 1; jcol <= k - 1; jcol++ )
      {
        if ( irow == iarray[jcol-1] + k - jcol || 
             irow == iarray[jcol-1] - ( k - jcol ) )
        {
          diag = true;
        }
      }

      if ( !diag )
      {
        ncan = ncan + 1;
        istack[nstack] = irow;
        nstack = nstack + 1;
      }
    }
  }

  istack[nstack] = ncan;
  nstack = nstack + 1;

  return;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
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
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( - x + 0.5 );
  }
  else
  {
    value =   ( int ) (  x + 0.5 );
  }
  return value;
}
//****************************************************************************80

float r4_uniform ( float a, float b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNIFORM returns a scaled pseudorandom R4.
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
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, float R4_UNIFORM, a number strictly between A and B.
//
{
  int i4_huge = 2147483647;
  int k;
  float value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }

  value = ( float ) ( *seed ) * 4.656612875E-10;

  value = a + ( b - a ) * value;

  return value;
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

double r8_add ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ADD adds two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the numbers to be added.
//
//    Output, double R8_ADD, the sum of X and Y.
//
{
  double value;

  value = x + y;

  return value;
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

double r8_gamma_log ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
//
//  Discussion:
//
//    The C MATH library includes a function LGAMMA ( X ) which should be
//    invoked instead of this function.
//
//    The program uses rational functions that theoretically approximate
//    LOG(GAMMA(X)) to at least 18 significant decimal digits.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by William Cody, Laura Stoltz.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    William Cody and Kenneth Hillstrom,
//    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
//    Mathematics of Computation,
//    Volume 21, 1967, pages 198-203.
//
//    Kenneth Hillstrom,
//    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
//    May 1969.
//
//    Hart, Et. Al.,
//    Computer Approximations,
//    Wiley and sons, New York, 1968.
//
//  Parameters:
//
//    Input, double X, the argument of the Gamma function.  X must be positive.
//
//    Output, double R8_GAMMA_LOG, the logarithm of the Gamma function of X.
//    If X <= 0.0, or if overflow would occur, the program returns the
//    value XINF, the largest representable floating point number.
//
//  Local Parameters:
//
//  BETA   - radix for the floating-point representation.
//
//  MAXEXP - the smallest positive power of BETA that overflows.
//
//  XBIG   - largest argument for which LN(GAMMA(X)) is representable
//           in the machine, i.e., the solution to the equation
//             LN(GAMMA(XBIG)) = BETA**MAXEXP.
//
//  FRTBIG - Rough estimate of the fourth root of XBIG
//
//
//  Approximate values for some important machines are:
//
//                            BETA      MAXEXP         XBIG
//
//  CRAY-1        (S.P.)        2        8191       9.62E+2461
//  Cyber 180/855
//    under NOS   (S.P.)        2        1070       1.72E+319
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)        2         128       4.08E+36
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)        2        1024       2.55D+305
//  IBM 3033      (D.P.)       16          63       4.29D+73
//  VAX D-Format  (D.P.)        2         127       2.05D+36
//  VAX G-Format  (D.P.)        2        1023       1.28D+305
//
//
//                           FRTBIG
//
//  CRAY-1        (S.P.)   3.13E+615
//  Cyber 180/855
//    under NOS   (S.P.)   6.44E+79
//  IEEE (IBM/XT,
//    SUN, etc.)  (S.P.)   1.42E+9
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)   2.25D+76
//  IBM 3033      (D.P.)   2.56D+18
//  VAX D-Format  (D.P.)   1.20D+9
//  VAX G-Format  (D.P.)   1.89D+76
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
  double corr;
  double d1 = - 5.772156649015328605195174E-01;
  double d2 =   4.227843350984671393993777E-01;
  double d4 =   1.791759469228055000094023;
  double frtbig = 1.42E+09;
  int i;
  double p1[8] = {
    4.945235359296727046734888, 
    2.018112620856775083915565E+02, 
    2.290838373831346393026739E+03, 
    1.131967205903380828685045E+04, 
    2.855724635671635335736389E+04, 
    3.848496228443793359990269E+04, 
    2.637748787624195437963534E+04, 
    7.225813979700288197698961E+03 };
  double p2[8] = {
    4.974607845568932035012064, 
    5.424138599891070494101986E+02, 
    1.550693864978364947665077E+04, 
    1.847932904445632425417223E+05, 
    1.088204769468828767498470E+06, 
    3.338152967987029735917223E+06, 
    5.106661678927352456275255E+06, 
    3.074109054850539556250927E+06 };
  double p4[8] = {
    1.474502166059939948905062E+04, 
    2.426813369486704502836312E+06, 
    1.214755574045093227939592E+08, 
    2.663432449630976949898078E+09, 
    2.940378956634553899906876E+010,
    1.702665737765398868392998E+011,
    4.926125793377430887588120E+011, 
    5.606251856223951465078242E+011 };
  double pnt68 = 0.6796875;
  double q1[8] = {
    6.748212550303777196073036E+01, 
    1.113332393857199323513008E+03, 
    7.738757056935398733233834E+03, 
    2.763987074403340708898585E+04, 
    5.499310206226157329794414E+04, 
    6.161122180066002127833352E+04, 
    3.635127591501940507276287E+04, 
    8.785536302431013170870835E+03 };
  double q2[8] = {
    1.830328399370592604055942E+02, 
    7.765049321445005871323047E+03, 
    1.331903827966074194402448E+05, 
    1.136705821321969608938755E+06, 
    5.267964117437946917577538E+06, 
    1.346701454311101692290052E+07, 
    1.782736530353274213975932E+07, 
    9.533095591844353613395747E+06 };
  double q4[8] = {
    2.690530175870899333379843E+03, 
    6.393885654300092398984238E+05, 
    4.135599930241388052042842E+07, 
    1.120872109616147941376570E+09, 
    1.488613728678813811542398E+010, 
    1.016803586272438228077304E+011, 
    3.417476345507377132798597E+011, 
    4.463158187419713286462081E+011 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double xbig = 4.08E+36;
  double xden;
  double xm1;
  double xm2;
  double xm4;
  double xnum;
  double xsq;
//
//  Return immediately if the argument is out of range.
//
  if ( x <= 0.0 || xbig < x )
  {
    return r8_huge ( );
  }

  if ( x <= r8_epsilon ( ) )
  {
    res = - log ( x );
  }
  else if ( x <= 1.5 )
  {
    if ( x < pnt68 )
    {
      corr = - log ( x );
      xm1 = x;
    }
    else
    {
      corr = 0.0;
      xm1 = ( x - 0.5 ) - 0.5;
    }

    if ( x <= 0.5 || pnt68 <= x )
    {
      xden = 1.0;
      xnum = 0.0;

      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm1 + p1[i];
        xden = xden * xm1 + q1[i];
      }

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) );
    }
    else
    {
      xm2 = ( x - 0.5 ) - 0.5;
      xden = 1.0;
      xnum = 0.0;
      for ( i = 0; i < 8; i++ )
      {
        xnum = xnum * xm2 + p2[i];
        xden = xden * xm2 + q2[i];
      }

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) );

    }
  }
  else if ( x <= 4.0 )
  {
    xm2 = x - 2.0;
    xden = 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm2 + p2[i];
      xden = xden * xm2 + q2[i];
    }

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) );
  }
  else if ( x <= 12.0 )
  {
    xm4 = x - 4.0;
    xden = - 1.0;
    xnum = 0.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = xnum * xm4 + p4[i];
      xden = xden * xm4 + q4[i];
    }

    res = d4 + xm4 * ( xnum / xden );
  }
  else
  {
    res = 0.0;

    if ( x <= frtbig )
    {

      res = c[6];
      xsq = x * x;

      for ( i = 0; i < 6; i++ )
      {
        res = res / xsq + c[i];
      }

    }

    res = res / x;
    corr = log ( x );
    res = res + sqrtpi - 0.5 * corr;
    res = res + x * ( corr - 1.0 );

  }

  return res;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r8_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r8_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

void r8vec_backtrack ( int n, int maxstack, double stack[], double x[], 
  int *indx, int *k, int *nstack, int ncan[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BACKTRACK supervises a backtrack search for a real vector.
//
//  Discussion:
//
//    The routine tries to construct a real vector one index at a time,
//    using possible candidates as supplied by the user.
//
//    At any time, the partially constructed vector may be discovered to be
//    unsatisfactory, but the routine records information about where the
//    last arbitrary choice was made, so that the search can be
//    carried out efficiently, rather than starting out all over again.
//
//    First, call the routine with INDX = 0 so it can initialize itself.
//
//    Now, on each return from the routine, if INDX is:
//      1, you've just been handed a complete candidate vector;
//         Admire it, analyze it, do what you like.
//      2, please determine suitable candidates for position X(K).
//         Return the number of candidates in NCAN(K), adding each
//         candidate to the end of STACK, and increasing NSTACK.
//      3, you're done.  Stop calling the routine;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2004
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
//    Input, int N, the number of positions to be filled in the vector.
//
//    Input, int MAXSTACK, the maximum length of the stack.
//
//    Input, double STACK[MAXSTACK], a list of all current candidates for
//    all positions 1 through K.
//
//    Input/output, double X[N], the partial or complete candidate vector.
//
//    Input/output, int *INDX, a communication flag.
//    On input,
//      0 to start a search.
//    On output:
//      1, a complete output vector has been determined and returned in X(1:N);
//      2, candidates are needed for position X(K);
//      3, no more possible vectors exist.
//
//    Inout/output, int *K, if INDX=2, the current vector index being considered.
//
//    Input/output, int *NSTACK, the current length of the stack.
//
//    Input/output, int NCAN[N], lists the current number of candidates for
//    positions 1 through K.
//
{
//
//  If this is the first call, request a candidate for position 1.
//
  if ( *indx == 0 )
  {
    *k = 1;
    *nstack = 0;
    *indx = 2;
    return;
  }
//
//  Examine the stack.
//
  for ( ; ; )
  {
//
//  If there are candidates for position K, take the first available
//  one off the stack, and increment K.
//
//  This may cause K to reach the desired value of N, in which case
//  we need to signal the user that a complete set of candidates
//  is being returned.
//
    if ( 0 < ncan[(*k)-1] )
    {
      x[(*k)-1] = stack[(*nstack)-1];
      *nstack = *nstack - 1;

      ncan[(*k)-1] = ncan[(*k)-1] - 1;

      if ( *k != n )
      {
        *k = *k + 1;
        *indx = 2;
      }
      else
      {
        *indx = 1;
      }

      break;
    }
//
//  If there are no candidates for position K, then decrement K.
//  If K is still positive, repeat the examination of the stack.
//
    else
    {
      *k = *k - 1;

      if ( *k <= 0 )
      {
        *indx = 3;
        break;
      }

    }

  }

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

int rgf_check ( int m, int f[] )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_CHECK checks a restricted growth function.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.
// 
//    Input, int F(M), the restricted growth function.
// 
//    Output, int IERROR, error flag.
//    0, no error.
//    -1, M is illegal.
//    I, entry I of the restricted growth function is illegal.
// 
{
  int fmax;
  int i;
  int ierror;

  ierror = 0;

  if ( m <= 0 )
  {
    ierror = -1;
    return ierror;
  }

  fmax = 0;
  for ( i = 0; i < m; i++ )
  {
    if ( f[i] <= 0 || fmax + 1 < f[i] )
    {
      ierror = i + 1;
      return ierror;
    }
    fmax = i4_max ( fmax, f[i] );
  }

  return ierror;
}
//****************************************************************************80

int rgf_enum ( int m )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_ENUM enumerates the restricted growth functions on M.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.  However, for the enumeration routine
//    only, it is legal to call with any value of M.
// 
//    Output, int RGF_ENUM, the number of restricted growth
//    functions.
// 
{
  int *b;
  int i;
  int j;
  int value;

  if ( m < 0 )
  {
    value = 0;
  }
  else if ( m == 0 )
  {
    value = 1;
  }
  else
  {
    b = new int[m+1];
    b[0] = 1;
    for ( j = 1; j <= m; j++ )
    {
      b[j] = 0;
      for ( i = 0; i < j; i++ )
      {
        b[j] = b[j] + i4_choose ( j - 1, i ) * b[i];
      }
    }
    value = b[m];

    delete [] b;
  }
  return value;
}
//****************************************************************************80

int *rgf_g_table ( int m )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_G_TABLE tabulates the generalized restricted growth functions.
// 
//  Example:
// 
//    M = 6
// 
//    D =  1    1    1    1    1    1    1
//         1    2    3    4    5    6    0
//         2    5   10   17   26    0    0
//         5   15   37   77    0    0    0
//        15   52  151    0    0    0    0
//        52  203    0    0    0    0    0
//       203    0    0    0    0    0    0
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int M, indicates how many rows and columns are to
//    be computed.  M must be nonnegative.
// 
//    Output, int RGF_G_TABLE[(M+1)*(M+1)], the first M+1 rows and
//    M+1 columns of the table of the number of generalized restricted growth
//    functions.  D(I,J) is the number of GRGF''s of length I with restriction
//    parameter J.
// 
{
  int *d;
  int i;
  int j;

  d = new int[(m+1)*(m+1)];

  for ( j = 0; j <= m; j++ )
  {
    d[0+j*(m+1)] = 1;
  }

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 0; j <= m; j++ )
    {
      if ( j <= m - i )
      {
        d[i+j*(m+1)] = j * d[i-1+j*(m+1)] + d[i-1+(j+1)*(m+1)];
      }
      else
      {
        d[i+j*(m+1)] = 0;
      }
    }
  }
  return d;
}
//****************************************************************************80

int rgf_rank ( int m, int f[] )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_RANK ranks a restricted growth function.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.
// 
//    Input, int F[M], the restricted growth function.
// 
//    Output, int RGF_RANK, the rank of the restricted growth
//    function.
// 
{
  int *d;
  int i;
  int ierror;
  int j;
  int rank;
// 
//  Check.
// 
  ierror = rgf_check ( m, f );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "RGF_RANK - Fatal error!\n";
    cout << "  The input array is illegal!\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Get the generalized restricted growth function table.
// 
  d = rgf_g_table ( m );

  rank = 0;
  j = 1;
  for ( i = 2; i <= m; i++ )
  {
    rank = rank + ( f[i-1] - 1 ) * d[m-i+j*(m+1)];
    j = i4_max ( j, f[i-1] );
  }

  delete [] d;

  return rank;
}
//****************************************************************************80

void rgf_successor ( int m, int f[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_SUCCESSOR generates the next restricted growth function.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.
// 
//    Input/output, int F[M], the restricted growth function.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
// 
{
  int fmax;
  int i;
  int ierror;
  int j;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    for ( i = 0; i < m; i++ )
    { 
       f[i] = 1;
    }
    rank = 0;
    return;
  }
// 
//  Check.
// 
  ierror = rgf_check ( m, f );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "RGF_SUCCESSOR - Fatal error!\n";
    cout << "  The input array is illegal!\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Find the first position from the right which can be incremented.
// 
  for ( i = m; 2 <= i; i-- )
  {
    fmax = 1;
    for ( j = 2; j < i; j++ )
    {
      fmax = i4_max ( fmax, f[j-1] );
    }
// 
//  Increment the function at this position, and set later entries to 1.
// 
    if ( f[i-1] != fmax + 1 )
    {
      f[i-1] = f[i-1] + 1;
      for ( j = i + 1; j <= m; j++ )
      {
        f[j-1] = 1;
      }
      rank = rank + 1;
      return;
    }
  }
// 
//  The final element was input.
//  Return the first element.
// 
  for ( i = 0; i < m; i++ )
  {
    f[i] = 1;
  }
  rank = 0;

  return;
}
//****************************************************************************80

void rgf_to_setpart ( int m, int f[], int &nsub, int s[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_TO_SETPART converts a restricted growth function to a set partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.
// 
//    Input, int F[M], the restricted growth function.
// 
//    Output, int NSUB, the number of nonempty subsets into
//    which the set is partitioned.
// 
//    Output, int S[M], describes the partition of a set of
//    M objects into NSUB nonempty subsets.  If element I of the
//    superset belongs to subset J, then S(I) = J.
// 
//    Output, int INDEX[M], lists the location in S of the last
//    element of each subset.  Thus, the elements of subset 1
//    are S(1) through S(INDEX(1)), the elements of subset 2
//    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
// 
{
  int i;
  int ierror;
  int j;
  int k;
// 
//  Check.
// 
  ierror = rgf_check ( m, f );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "RGF_TO_SETPART - Fatal error!\n";
    cout << "  The input array is illegal!\n";
    cout << "  IERROR = " << ierror << "\n";
    exit ( 1 );
  }
// 
//  Determine the number of subsets.
// 
  nsub = i4vec_max ( m, f );
// 
//  Initialize.
//
  for ( i = 0; i < m; i++ )
  {
    s[i] = 0;
  }
  for ( i = 0; i < nsub; i++ )
  {
    index[i] = 0;
  }
// 
//  For each subset I, collect the indices of F which have value I.
//  These are the elements of the I-th subset.
// 
  k = 0;
  for ( i = 1; i <= nsub; i++ )
  {
    for ( j = 1; j <= m; j++ )
    {
      if ( f[j-1] == i )
      {
        k = k + 1;
        s[k-1] = j;
      }
    }
    index[i-1] = k;
  }
  return;
}
//****************************************************************************80

int *rgf_unrank ( int rank, int m )

//****************************************************************************80
// 
//  Purpose:
//
//    RGF_UNRANK returns the restricted growth function of a given rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the restricted growth
//    function.
// 
//    Input, int M, the domain of the RGF is the integers
//    from 1 to M.  M must be positive.
// 
//    Output, int RGF_UNRANK[M], the restricted growth function.
// 
{
  int *d;
  int *f;
  int i;
  int j;
  int nrgf;
  int rank_copy;
// 
//  Check.
// 
  if ( m < 1 )
  {
    cout << "\n";
    cout << "RGF_UNRANK - Fatal error!\n";
    cout << "  Input M is illegal.\n";
    exit ( 1 );
  }

  nrgf = rgf_enum ( m );

  if ( rank < 0 || nrgf < rank )
  {
    cout << "\n";
    cout << "RGF_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }
// 
//  Get the generalized restricted growth function table.
// 
  d = rgf_g_table ( m );

  f = new int[m];

  rank_copy = rank;
  j = 1;
  f[0] = 1;

  for ( i = 2; i <= m; i++ )
  {
    if ( j * d[m-i+j*(m+1)] <= rank_copy )
    {
      f[i-1] = j + 1;
      rank_copy = rank_copy - j * d[m-i+j*(m+1)];
      j = j + 1;
    }
    else
    {
      f[i-1] = 1 + ( rank_copy / d[m-i+j*(m+1)] );
      rank_copy = rank_copy % d[m-i+j*(m+1)];
    }
  }

  delete [] d;

  return f;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int setpart_check ( int m, int nsub, int s[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SETPART_CHECK checks a set partition.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 August 2012
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
//    Input, int M, the number of elements of the set.
//    M must be positive.
// 
//    Input, int NSUB, the number of nonempty subsets into
//    which the set is partitioned.  1 <= NSUB <= M.
//
//    Input, int S[M], contains the integers from 1 to M,
//    grouped into subsets as described by INDEX.
// 
//    Input, int INDEX[NSUB], lists the location in S of the
//    last element of each subset.  Thus, the elements of subset 1
//    are S(1) through S(INDEX(1)), the elements of subset 2
//    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
// 
//    Output, int SETPART_CHECK, error flag.
//    0, no error.
//    -I, the I-th element of INDEX is illegal.
//    +I, the I-th element of S is illegal.
// 
{
  int i;
  int ierror;
  int imin;
  int j;

  ierror = 0;
// 
//  Check INDEX.
// 
  imin = 0;
  for ( i = 0; i < nsub; i++ )
  {
    if ( index[i] <= imin || m < index[i] )
    {
      ierror = - ( i + 1 );
      return ierror;
    }
    imin = index[i];
  }
// 
//  Check the elements of S.
// 
  for ( i = 0; i < nsub; i++ )
  {
    if ( s[i] <= 0 || m < s[i] )
    {
      ierror = i + 1;
      return ierror;
    }

    for ( j = 0; j < i; j++ )
    {
      if ( s[j] == s[i] )
      {
        ierror = i + 1;
        return ierror;
      }
    }
  }

  return ierror;
}
//****************************************************************************80

int setpart_enum ( int m )

//****************************************************************************80
// 
//  Purpose:
//
//    SETPART_ENUM enumerates the partitions of a set of M elements.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    27 July 2011
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
//    Input, int M, the number of elements in the set.
//    M must be positive.  However, for the enumeration routine only,
//    it is legal to call with any value of M.
// 
//    Output, int SETPART_ENUM, the number of partitions of the set.
// 
{
  int *b;
  int i;
  int j;
  int value;

  if ( m < 0 )
  {
    value = 0;
  }
  else if ( m == 0 )
  {
    value = 1;
  }
  else
  {
    b = new int[m+1];
    b[0] = 1;
    for ( j = 1; j <= m; j++ )
    {
      b[j] = 0;
      for ( i = 0; i < j; i++ )
      {
        b[j] = b[j] + i4_choose ( j - 1, i ) * b[i];
      }
    }
    value = b[m];

    delete [] b;
  }

  return value;
}
//****************************************************************************80

int *setpart_to_rgf ( int m, int nsub, int s[], int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SETPART_TO_RGF converts a set partition to a restricted growth function.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int M, the number of elements of the set.
//    M must be positive.
// 
//    Input, int NSUB, the number of nonempty subsets into
//    which the set is partitioned.  1 <= NSUB <= M.
//
//    Input, int S(M), contains the integers from 1 to M,
//    grouped into subsets as described by INDEX.
// 
//    Input, int INDEX(NSUB), lists the location in S of the
//    last element of each subset.  Thus, the elements of subset 1
//    are S(1) through S(INDEX(1)), the elements of subset 2
//    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
// 
//    Output, int SETPART_TO_RGF[M], the restricted growth function from
//    M to NSUB.
// 
{
  int *f;
  int i;
  int ierror;
  int k;
  int khi;
  int klo;
// 
//  Check.
// 
  ierror = setpart_check ( m, nsub, s, index );

  if ( ierror != 0 )
  {
    cout << "\n";
    cout << "SETPART_TO_RGF - Fatal error!\n";
    cout << "  The input array is illegal.\n";
    exit ( 1 );
  }

  f = new int[m];

  khi = 0;
  for ( i = 1; i <= nsub; i++ )
  {
    klo = khi + 1;
    khi = index[i-1];
    for ( k = klo; k <= khi; k++ )
    {
      f[s[k-1]-1] = i;
    }
  }
  return f;
}
//****************************************************************************80

int *stirling_numbers1 ( int m, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    STIRLING_NUMBERS1 computes Stirling numbers of the first kind.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int M, the maximum row to compute.
//    M must be nonnegative.
// 
//    Input, int N, the maximum column to compute.
//    N must be nonnegative.
// 
//    Output, int S[(M+1)*(N+1)], the first M+1 rows and N+1 columns
//    of the table of Stirling numbers of the first kind.
// 
{
  int i;
  int j;
  int *s;

  s = new int[(m+1)*(n+1)];

  s[0+0*(m+1)] = 1;
  for ( j = 1; j <= n; j++ )
  {
    s[0+j*(m+1)] = 0;
  }

  for ( i = 1; i <= m; i++ )
  {
    s[i+0*(m+1)] = 0;
  }
// 
//  This loop may be extraneous.
// 
  for ( i = 0; i <= i4_min ( m, n - 1 ); i++ )
  {
    s[i+(i+1)*(m+1)] = 0;
  }

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j <= i )
      {
        s[i+j*(m+1)] = s[i-1+(j-1)*(m+1)] - ( i - 1 ) * s[i-1+j*(m+1)];
      }
      else
      {
        s[i+j*(m+1)] = 0;
      }
    }
  }

  return s;
}
//****************************************************************************80

int *stirling_numbers2 ( int m, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    STIRLING_NUMBERS2 computes Stirling numbers of the second kind.
// 
//  Discussion:
// 
//    The reference has a typographical error, referring to
//    S(I-J,J-1) instead of S(I-1,J-1).
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int M, the maximum row to compute.
//    M must be nonnegative.
// 
//    Input, int N, the maximum column to compute.
//    N must be nonnegative.
// 
//    Output, int S[(M+1)*(N+1)], the first M+1 rows and N+1 columns
//    of the table of Stirling numbers of the second kind.
// 
{
  int i;
  int j;
  int *s;

  s = new int[(m+1)*(n+1)];

  s[0+0*(m+1)] = 1;
  for ( j = 1; j <= n; j++ )
  {
    s[0+j*(m+1)] = 0;
  }
  for ( i = 1; i <= m; i++ )
  {
    s[i+0*(m+1)] = 0;
  }
// 
//  This loop may be extraneous.
// 
  for ( i = 0; i <= i4_min ( m, n - 1 ); i++ )
  {
    s[i+(i+1)*(m+1)] = 0;
  }

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      if ( j <= i )
      {
        s[i+j*(m+1)] = j * s[i-1+j*(m+1)] + s[i-1+(j-1)*(m+1)];
      }
      else
      {
        s[i+j*(m+1)] = 0;
      }
    }
  }

  return s;
}
//****************************************************************************80

void subset_check ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_CHECK checks a subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    28 July 2011
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
//    Input, int T[N], the subset.  If T(I) = 0, item I is
//    not in the subset; if T(I) = 1, item I is in the subset.
// 
{
  int i;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "SUBSET_CHECK - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] != 0 && t[i] != 1 )
    {
      cerr << "\n";
      cerr << "SUBSET_CHECK - Fatal error!\n";
      cerr << "  T[I] = 0 or 1 is required.\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************80

int subset_colex_rank ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_COLEX_RANK computes the colexicographic rank of a subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    22 August 2011
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
//    Input, int N, the number of items in the master set.
//    N must be positive.
// 
//    Input, int T[N], the subset.  If T(I) = 0, item I is
//    not in the subset; if T(I) = 1, item I is in the subset.
// 
//    Output, int SUBSET_COLEX_RANK, the rank of the subset.
// 
{
  int i;
  int rank;
// 
//  Check.
// 
  subset_check ( n, t );

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] == 1 )
    {
      rank = rank + i4_power ( 2, i );
    }
  }
  return rank;
}
//****************************************************************************80

void subset_colex_successor ( int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_COLEX_SUCCESSOR computes the subset colexicographic successor.
// 
//  Discussion:
// 
//    In the original code, there is a last element with no successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    If the input T was the last in the ordering, then the output T
//    will be the first.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
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
// 
//  Check.
// 
  subset_check ( n, t );

  for ( i = 0; i < n; i++ )
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
  rank = 0;

  return;
}
//****************************************************************************80

int *subset_colex_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_COLEX_UNRANK computes the subset of given colexicographic rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the subset.
// 
//    Input, int N, the number of items in the master set.
//    N must be positive.
// 
//    Output, int SUBSET_COLEX_UNRANK[N], the subsetof the given rank.
//    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is
//    in the subset.
// 
{
  int i;
  int nsub;
  int rank_copy;
  int *t;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "SUBSET_COLEX_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  nsub = subset_enum ( n );

  if ( rank < 0 || nsub < rank )
  {
    cout << "\n";
    cout << "SUBSET_COLEX_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  rank_copy = rank;
  t = new int[n];

  for ( i = 0; i < n; i++ )
  {
    if ( ( rank_copy % 2 ) == 1 )
    {
      t[i] = 1;
    }
    else
    {
      t[i] = 0;
    }
    rank_copy = rank_copy / 2;
  }
  return t;
}
//****************************************************************************80

int *subset_complement ( int n, int a[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_COMPLEMENT computes the complement of a set.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the order of the master set, of which A is
//    a subset.  N must be positive.
// 
//    Input, int A[N], a subset of the master set.
//    A(I) = 0 if the I-th element is in the subset A, and is
//    1 otherwise.
// 
//    Output, int SUBSET_COMPLEMENT[N], the complement of A.
// 
{
  int *b;
  int i;
// 
//  Check.
// 
  subset_check ( n, a );

  b = new int[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = 1 - a[i];
  }
  return b;
}
//****************************************************************************80

int subset_distance ( int n, int t1[], int t2[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_DISTANCE computes the Hamming distance between two sets.
// 
//  Discussion:
// 
//    The sets T1 and T2 are assumed to be subsets of a set of N elements.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the order of the master set, of which T1 and
//    T2 are subsets.  N must be positive.
// 
//    Input, int T1[N], T2[N], two subsets of the master set.
//    T1(I) = 0 if the I-th element is in the subset T1, and is
//    1 otherwise; T2 is defined similarly.
// 
//    Output, int SUBSET_DISTANCE, the Hamming distance between T1 and T2,
//    defined as the number of elements of the master set which are
//    in either T1 or T2 but not both.
// 
{
  int dist;
  int i;
// 
//  Check.
// 
  subset_check ( n, t1 );

  subset_check ( n, t2 );

  dist = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( ( t1[i] == 0 && t2[i] != 0 ) || ( t1[i] != 0 && t2[i] == 0 ) )
    {
      dist = dist + 1;
    }
  }
  return dist;
}
//****************************************************************************80

int subset_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_ENUM enumerates the subsets of a set with N elements.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of elements in the set.
//    N must be at least 0.
// 
//    Output, int SUBSET_ENUM, the number of distinct elements.
// 
{
  int value;

  value = i4_power ( 2, n );

  return value;
}
//****************************************************************************80

int *subset_intersect ( int n, int a[], int b[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_INTERSECT computes the intersection of two sets.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the order of the master set, of which A and
//    B are subsets.  N must be positive.
// 
//    Input, int A[N], B[N], two subsets of the master set.
//    A(I) = 0 if the I-th element is in the subset A, and is
//    1 otherwise; B is defined similarly.
// 
//    Output, int SUBSET_INTERSECT[N], the intersection of A and B.
// 
{
  int *c;
  int i;
// 
//  Check.
// 
  subset_check ( n, a );

  subset_check ( n, b );

  c = new int[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = i4_min ( a[i], b[i] );
  }
  return c;
}
//****************************************************************************80

int subset_lex_rank ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_LEX_RANK computes the lexicographic rank of a subset.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    22 August 2011
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
//    Input, int N, the number of items in the master set.
//    N must be positive.
// 
//    Input, int T[N], the subset.  If T(I) = 0, item I is
//    not in the subset; if T(I) = 1, item I is in the subset.
// 
//    Output, int SUBSET_LEX_RANK, the rank of the subset.
// 
{
  int i;
  int rank;
// 
//  Check.
// 
  subset_check ( n, t );

  rank = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( t[i] == 1 )
    {
      rank = rank + i4_power ( 2, n - i - 1 );
    }
  }

  return rank;
}
//****************************************************************************80

void subset_lex_successor ( int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_LEX_SUCCESSOR computes the subset lexicographic successor.
// 
//  Discussion:
// 
//    In the original code, there is a last element with no successor.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    27 July 2011
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
//    If the input T was the last in the ordering, then the output T
//    will be the first.
// 
//    Input/output, int &RANK, the rank.
//    If RANK = -1 on input, then the routine understands that this is
//    the first call, and that the user wishes the routine to supply
//    the first element in the ordering, which has RANK = 0.
//    In general, the input value of RANK is increased by 1 for output,
//    unless the very last element of the ordering was input, in which
//    case the output value of RANK is 0.
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
// 
//  Check.
// 
  subset_check ( n, t );

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
  rank = 0;

  return;
}
//****************************************************************************80

int *subset_lex_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_LEX_UNRANK computes the subset of given lexicographic rank.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int RANK, the rank of the subset.
// 
//    Input, int N, the number of items in the master set.
//    N must be positive.
// 
//    Output, int SUBSET_LEX_UNRANK[N], the subset of the given rank.
//    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is in
//    the subset.
// 
{
  int i;
  int nsub;
  int rank_copy;
  int *t;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "SUBSET_LEX_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  nsub = subset_enum ( n );

  if ( rank < 0 || nsub < rank )
  {
    cout << "\n";
    cout << "SUBSET_LEX_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }

  rank_copy = rank;
  t = new int[n];

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( ( rank_copy % 2 ) == 1 )
    {
      t[i] = 1;
    }
    else
    {
      t[i] = 0;
    }
    rank_copy = rank_copy / 2;
  }

  return t;
}
//****************************************************************************80

int *subset_union ( int n, int a[], int b[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_UNION computes the union of two sets.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the order of the master set, of which A and
//    B are subsets.  N must be positive.
// 
//    Input, int A[N], B[N], two subsets of the master set.
//    A(I) = 0 if the I-th element is in the subset A, and is
//    1 otherwise; B is defined similarly.
// 
//    Output, int SUBSET_UNION[N], the union of A and B.
// 
{
  int *c;
  int i;
// 
//  Check.
// 
  subset_check ( n, a );

  subset_check ( n, b );

  c = new int[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = i4_max ( a[i], b[i] );
  }

  return c;
}
//****************************************************************************80

int subset_weight ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_WEIGHT computes the Hamming weight of a set.
// 
//  Discussion:
// 
//    The Hamming weight is simply the number of elements in the set.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the order of the master set, of which T
//    is a subset.  N must be positive.
// 
//    Input, int T[N], defines the subset T.
//    T(I) is 1 if I is an element of T, and 0 otherwise.
// 
//    Output, int SUBSET_WEIGHT, the Hamming weight of the subset T.
// 
{
  int weight;
// 
//  Check.
// 
  subset_check ( n, t );

  weight = i4vec_sum ( n, t );

  return weight;
}
//****************************************************************************80

int *subset_xor ( int n, int a[], int b[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSET_XOR computes the symmetric difference of two sets.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the order of the master set, of which A and
//    B are subsets.  N must be positive.
// 
//    Input, int A[N], B[N], two subsets of the master set.
//    A(I) = 0 if the I-th element is in the subset A, and is
//    1 otherwise; B is defined similarly.
// 
//    Output, int SUBSET_XOR[N], the symmetric difference of A and B.
// 
{
  int *c;
  int i;
// 
//  Check.
// 
  subset_check ( n, a );

  subset_check ( n, b );

  c = new int[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = i4_max ( a[i], b[i] ) - i4_min ( a[i], b[i] );
  }

  return c;
}
//****************************************************************************80

int subsetsum_swap ( int n, int a[], int sum_desired, int index[] )

//****************************************************************************80
// 
//  Purpose:
//
//    SUBSETSUM_SWAP seeks a solution of the subset sum problem by swapping.
// 
//  Discussion:
// 
//    Given a collection of N not necessarily distinct positive integers A(I),
//    and a positive integer SUM_DESIRED, select a subset of the values so that
//    their sum is as close as possible to SUM_DESIRED without exceeding it.
// 
//  Algorithm:
// 
//    Start with no values selected, and SUM_ACHIEVED = 0.
// 
//    Consider each element A(I):
// 
//      If A(I) is not selected and SUM_ACHIEVED + A(I) <= SUM_DESIRED,
//        select A(I).
// 
//      If A(I) is still not selected, and there is a selected A(J)
//      such that SUM_GOT < SUM_ACHIEVED + A(I) - A(J),
//        select A(I) and deselect A(J).
// 
//      If no items were selected on this sweep,
//        exit.
//      Otherwise,
//        repeat the search.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of values.  N must be positive.
// 
//    Input/output, int A[N], a collection of positive values.
//    On output, A has been sorted into descending order.
// 
//    Input, int SUM_DESIRED, the desired sum.
// 
//    Output, int INDEX[N]; INDEX(I) is 1 if A(I) is part of the
//    sum, and 0 otherwise.
// 
//    Output, int SUBSETSUM_SWAP, the sum of the selected
//    elements.
// 
{
  int i;
  int j;
  int nmove;
  int sum_achieved;
// 
//  Initialize.
// 
  sum_achieved = 0;

  for ( i = 0; i < n; i++ )
  {
    index[i] = 0;
  }
// 
//  Sort into descending order.
// 
  i4vec_sort_insert_d ( n, a );

  for ( ; ; )
  {
    nmove = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( index[i] == 0 )
      {
        if ( sum_achieved + a[i] <= sum_desired )
        {
          index[i] = 1;
          sum_achieved = sum_achieved + a[i];
          nmove = nmove + 1;
          continue;
        }
      }

      if ( index[i] == 0 )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( index[j] == 1 )
          {
            if ( sum_achieved < sum_achieved + a[i] - a[j] &&
              sum_achieved + a[i] - a[j] <= sum_desired )
            {
              index[j] = 0;
              index[i] = 1;
              nmove = nmove + 2;
              sum_achieved = sum_achieved + a[i] - a[j];
              break;
            }
          }
        }
      }
    }
    if ( nmove <= 0 )
    {
      break;
    }
  }

  return sum_achieved;
}
//****************************************************************************80

void tableau_check ( int n, int tab[] )

//****************************************************************************80
// 
//  Purpose:
//
//    TABLEAU_CHECK checks a 2 by N tableau.
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
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
// 
//  Parameters:
// 
//    Input, int N, the number of columns in the tableau.
//    N must be positive.
// 
//    Input, int TAB[2*N], a 2 by N tableau.
// 
{
  int i;
  int j;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "TABLEAU_CHECK - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
// 
//  The entries must be between 0 and 2*N.
// 
  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( tab[i+j*2] < 1 || 2 * n < tab[i+j*2] )
      {
        cerr << "\n";
        cerr << "TABLEAU_CHECK - Fatal error!\n";
        cerr << "  TAB[I][J] < 1 or N < TAB[I][J].\n";
        exit ( 1 );
      }
    }
  }
// 
//  The entries must be increasing to the right.
//
  for ( i = 0; i < 2; i++ )
  {
    for ( j = 1; j < n; j++ )
    {
      if ( tab[i+j*2] <= tab[i+(j-1)*2] )
      {
        cerr << "\n";
        cerr << "TABLEAU_CHECK - Fatal error!\n";
        cerr << "  TAB[I][J] < TAB[I][J-1].\n";
        exit ( 1 );
      }
    }
  }
// 
//  The entries must be increasing down.
// 
  i = 1;
  for ( j = 0; j < n; j++ )
  {
    if ( tab[i+j*2] <= tab[i-1+j*2] )
    {
      cerr << "\n";
      cerr << "TABLEAU_CHECK - Fatal error!\n";
      cerr << "  TAB[I][J] <= TAB[I-1][J].\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************80

int tableau_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    TABLEAU_ENUM enumerates the 2 by N standard tableaus.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
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
//    Input, int N, the number of columns in the tableau.
//    N must be nonnegative.
// 
//    Output, int TABLEAU_ENUM, the number of 2 by N standard tableaus.
// 
{
  int value;

  value = i4_choose ( 2 * n, n ) / ( n + 1 );

  return value;
}
//****************************************************************************80

int *tableau_to_bal_seq ( int n, int tab[] )

//****************************************************************************80
// 
//  Purpose:
//
//    TABLEAU_TO_BAL_SEQ converts a 2 by N tableau to a balanced sequence.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of 0's (and 1's) in the sequence.
//    N must be positive.
// 
//    Input, int TAB[2*N], a 2 by N tableau.
// 
//    Output, int TABLEAU_TO_BAL_SEQ[2*N], a balanced sequence.
// 
{
  int i;
  int j;
  int *t;
// 
//  Check.
// 
  tableau_check ( n, tab );

  t = new int[2*n];

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      t[tab[i+j*2]-1] = i;
    }
  }

  return t;
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
//****************************************************************************80

void tree_check ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_CHECK checks a tree.
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
//    Donald Kreher, Douglas Simpson,
//    Combinatorial Algorithms,
//    CRC Press, 1998,
//    ISBN: 0-8493-3988-X,
//    LC: QA164.K73.
// 
//  Parameters:
// 
//    Input, int N, the number of nodes in the tree.
//    N must be positive.
// 
//    Input, int T[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
{
  int *d;
  int i;
  int j;
  int k;
  int x;
  int y;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "TREE_CHECK - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n - 1; j++ )
    {
      if ( t[i+j*2] < 1 || n < t[i+j*2] )
      {
        cerr << "\n";
        cerr << "TREE_CHECK - Fatal error!\n";
        cerr << "  1 <= T[I][J] <= N is required.\n";
        exit ( 1 );
      }
    }
  }
// 
//  Compute the degree of each node.
//
  d = edge_degree ( n, n - 1, t );
// 
//  Delete a node of degree 1, N-1 times.
// 
  for ( k = 1; k <= n - 1; k++ )
  {
    x = 1;

    while ( d[x-1] != 1 )
    {
      x = x + 1;
      if ( n < x )
      {
        cerr << "\n";
        cerr << "TREE_CHECK - Fatal error!\n";
        cerr << "  Could not locate a node of degree 1.\n";
        exit ( 1 );
      }
    }
// 
//  Find its neighbor.
// 
    j = 1;

    for ( ; ; )
    {
      if ( t[0+(j-1)*2] == x )
      {
        y = t[1+(j-1)*2];
        break;
      }

      if ( t[1+(j-1)*2] == x )
      {
        y = t[0+(j-1)*2];
        break;
      }

      j = j + 1;

      if ( n < j )
      {
        cerr << "\n";
        cerr << "TREE_CHECK - Fatal error!\n";
        cerr << "  Hard to explain.\n";
        exit ( 1 );
      }
    }
// 
//  Delete the edge.
// 
    t[0+(j-1)*2] = - t[0+(j-1)*2];
    t[1+(j-1)*2] = - t[1+(j-1)*2];

    d[x-1] = d[x-1] - 1;
    d[y-1] = d[y-1] - 1;
  }

  for ( j = 0; j < n - 1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      t[i+j*2] = - t[i+j*2];
    }
  }

  delete [] d;

  return;
}
//****************************************************************************80

int tree_enum ( int n )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_ENUM enumerates the trees on N nodes.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    24 July 2011
// 
//  Author:
// 
//    John Burkardt
// 
//  Parameters:
// 
//    Input, int N, the number of nodes in each tree.
//    N must normally be at least 3, but for this routine,
//    any value of N is allowed.
// 
//    Output, int TREE_ENUM, the number of distinct elements.
// 
{
  int value;

  if ( n < 1 )
  {
    value = 0;
  }
  else if ( n == 1 )
  {
    value = 1;
  }
  else if ( n == 2 )
  {
    value = 1;
  }
  else
  {
    value = i4_power ( n, n - 2 );
  }
  return value;
}
//****************************************************************************80

int tree_rank ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_RANK ranks a tree.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input, int T[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
//    Output, int RANK, the rank of the tree.
// 
{
  int *p;
  int rank;
// 
//  Check the tree.
// 
  tree_check ( n, t );
// 
//  Convert the tree to a Pruefer code.
// 
  p = tree_to_pruefer ( n, t );
// 
//  Find the rank of the Pruefer code.
// 
  rank = pruefer_rank ( n, p );

  delete [] p;

  return rank;
}
//****************************************************************************80

void tree_successor ( int n, int t[], int &rank )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_SUCCESSOR returns the successor of a tree.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    26 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Input/output, int T[2*(N-1)], describes the edges of the
//    tree as pairs of nodes.  On output, the input tree has been replaced
//    by its successor.
// 
//    Input/output, int &RANK, the rank of the tree.
// 
{
  int i;
  int *p;
// 
//  Return the first element.
// 
  if ( rank == -1 )
  {
    p = new int[n-2];

    for ( i = 0; i < n - 2; i++ )
    {
      p[i] = 1;
    }
    pruefer_to_tree ( n, p, t );
    rank = 0;
    delete [] p;
    return;
  }
// 
//  Check the tree.
// 
  tree_check ( n, t );
// 
//  Convert the tree to a Pruefer code.
// 
  p = tree_to_pruefer ( n, t );
// 
//  Find the successor of the Pruefer code.
// 
  pruefer_successor ( n, p, rank );
// 
//  Convert the Pruefer code to the tree.
// 
  pruefer_to_tree ( n, p, t );

  delete [] p;

  return;
}
//****************************************************************************80

int *tree_to_pruefer ( int n, int t[] )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_TO_PRUEFER converts a tree to a Pruefer code.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    25 July 2011
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
//    Input, int N, the number of nodes in the tree.
//    N must be positive.
// 
//    Input, int T[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
//    Output, int TREE_TO_PRUEFER[N-2], the Pruefer code for the tree.
// 
{
  int *d;
  int i;
  int j;
  int k;
  int *p;
  int x;
  int y;
// 
//  Check.
// 
  tree_check ( n, t );
// 
//  Compute the degree of each node.
// 
  d = edge_degree ( n, n - 1, t );

  p = new int[n-2];

  for ( j = 1; j <= n - 2; j++ )
  {
// 
//  Find a node of degree 1.
// 
    x = n;
    while ( d[x-1] != 1 )
    {
      x = x - 1;
    }
// 
//  Find its neighbor.
// 
    k = 1;

    for ( ; ; )
    {
      if ( t[0+(k-1)*2] == x )
      {
        y = t[1+(k-1)*2];
        break;
      }

      if ( t[1+(k-1)*2] == x )
      {
        y = t[0+(k-1)*2];
        break;
      }
      k = k + 1;
    }
// 
//  Store the neighbor.
// 
    p[j-1] = y;
// 
//  Delete the edge from the tree.
// 
    d[x-1] = d[x-1] - 1;
    d[y-1] = d[y-1] - 1;

    t[0+(k-1)*2] = - t[0+(k-1)*2];
    t[1+(k-1)*2] = - t[1+(k-1)*2];
  }
// 
//  Remove the negative signs from the first N-2 columns of the tree.
// 
  for ( j = 0; j < n - 2; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      t[i+j*2] = - t[i+j*2];
    }
  }

  delete [] d;

  return p;
}
//****************************************************************************80

int *tree_unrank ( int rank, int n )

//****************************************************************************80
// 
//  Purpose:
//
//    TREE_UNRANK unranks a tree.
// 
//  Licensing:
// 
//    This code is distributed under the GNU LGPL license.
// 
//  Modified:
// 
//    27 July 2011
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
//    Input, int RANK, the rank of the tree.
// 
//    Input, int N, the number of nodes in the tree.
//    N must be at least 3.
// 
//    Output, int T[2*(N-1)], describes the edges of the tree
//    as pairs of nodes.
// 
{
  int *p;
  int *t;
  int tree_num;
// 
//  Check.
// 
  if ( n < 1 )
  {
    cout << "\n";
    cout << "TREE_UNRANK - Fatal error!\n";
    cout << "  Input N is illegal.\n";
    exit ( 1 );
  }

  tree_num = tree_enum ( n );

  if ( rank < 0 || tree_num < rank )
  {
    cout << "\n";
    cout << "TREE_UNRANK - Fatal error!\n";
    cout << "  The input rank is illegal.\n";
    exit ( 1 );
  }
// 
//  Unrank the Pruefer code.
// 
  p = pruefer_unrank ( rank, n );
// 
//  Convert the Pruefer code to a tree.
// 
  t = pruefer_to_tree_new ( n, p );

  delete [] p;

  return t;
}
