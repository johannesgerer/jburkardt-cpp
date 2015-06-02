# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "knapsack_01.hpp"

//****************************************************************************80

int *knapsack_01 ( int n, int w[], int c )

//****************************************************************************80
//
//  Purpose:
//
//    KNAPSACK_01 seeks a solution of the 0/1 Knapsack problem.
//
//  Discussion:
//
//    In the 0/1 knapsack problem, a knapsack of capacity C is given,
//    as well as N items, with the I-th item of weight W(I).
//
//    A selection is "acceptable" if the total weight is no greater than C.
//
//    It is desired to find an optimal acceptable selection, that is,
//    an acceptable selection such that there is no acceptable selection
//    of greater weight.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of weights.
//
//    Input, inte W[N], the weights.
//
//    Input, int C, the maximum weight.
//
//    Output, int KNAPSACK_01[N], is a binary vector which defines an 
//    optimal selection.  It is 1 for the weights to be selected, and 
//    0 otherwise.
//
{
  int i;
  int iadd;
  bool more;
  int ncard;
  int *s;
  int *s_test;
  int t;
  int t_test;

  s = new int[n];
  s_test = new int[n];

  more = false;
  ncard = 0;

  for ( i = 0; i < n; i++ )
  {
    s_test[i] = 0;
  }
  t_test = 0;

  for ( i = 0; i < n; i++ )
  {
    s[i] = s_test[i];
  }
  t = 0;

  for ( ; ; )
  {
    subset_gray_next ( n, s_test, more, ncard, iadd );
    t_test = 0;
    for ( i = 0; i < n; i++ )
    {
      t_test = t_test + s_test[i] * w[i];
    }

    if ( t < t_test && t_test <= c )
    {
      t = t_test;
      for ( i = 0; i < n; i++ )
      {
        s[i] = s_test[i];
      }
    }

    if ( ! more )
    {
      break;
    }
  }

  delete [] s_test;

  return s;
}
//****************************************************************************80

void subset_gray_next ( int n, int a[], bool &more, int &ncard, int &iadd )

//****************************************************************************80
//
//  Purpose:
//
//    SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
//
//  Discussion:
//
//    It generates the subsets one at a time, by adding or subtracting
//    exactly one element on each step.
//
//    The user should set MORE = .FALSE. and the value of N before
//    the first call.  On return, the user may examine A which contains
//    the definition of the new subset, and must check .MORE., because
//    as soon as it is .FALSE. on return, all the subsets have been
//    generated and the user probably should cease calling.
//
//    The first set returned is the empty set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 May 2003
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
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input/output, int A[N].  On each return, the Gray code for the newly
//    generated subset.  A[I] = 0 if element I is in the subset, 1 otherwise.
//
//    Input/output, bool &MORE.  Set this variable FALSE before
//    the first call.  Normally, MORE will be returned TRUE but once
//    all the subsets have been generated, MORE will be
//    reset FALSE on return and you should stop calling the program.
//
//    Input/output, int &NCARD, the cardinality of the set returned,
//    which may be any value between 0 (the empty set) and N (the
//    whole set).
//
//    Output, int &IADD, the element which was added or removed to the
//    previous subset to generate the current one.  Exception:
//    the empty set is returned on the first call, and IADD is set to -1.
{
  int i;
//
//  First set returned is the empty set.
//
  if ( !more )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    iadd = 0;
    ncard = 0;
    more = true;
  }
  else
  {
    iadd = 1;

    if ( ( ncard % 2 ) != 0 )
    {
      for ( ; ; )
      {
        iadd = iadd + 1;
        if ( a[iadd-2] != 0 )
        {
          break;
        }
      }
    }

    a[iadd-1] = 1 - a[iadd-1];
    ncard = ncard + 2 * a[iadd-1] - 1;
//
//  Last set returned is the singleton A(N).
//
    if ( ncard == a[n-1] )
    {
      more = false;
    }
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
