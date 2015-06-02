# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iostream>
# include <iomanip>

using namespace std;

# include "sort_rc.hpp"

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

void sort_rc ( int n, int &indx, int &i, int &j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_RC externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list of data is not passed to the routine.  Hence this
//    routine may be used to sort integers, reals, numbers, names,
//    dates, shoe sizes, and so on.  After each call, the routine asks
//    the user to compare or interchange two items, until a special
//    return value signals that the sorting is completed.
//
//    Note that this function uses internal persistent memory during the sort.
//
//  Example:
//
//    n = 100;
//    indx = 0;
//
//    for ( ; ; )
//    {
//      sort_rc ( n, indx, i, j, isgn );
//
//      if ( indx < 0 )
//      {
//        isgn = 1;
//        if ( a[i-1] <= a[j-1] )
//        {
//          isgn = -1;
//        }
//      }
//      else if ( 0 < indx )
//      {
//        k      = a[i-1];
//        a[i-1] = a[j-1];
//        a[j-1] = k;
//      }
//      else
//      {
//        break;
//      }
//    }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the length of the input list.
//
//    Input/output, int &INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int &I, &J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k_save = 0;
  static int l_save = 0;
  static int n_save = 0;
//
//  INDX = 0: This is the first call.
//
  if ( indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k_save = n / 2;
    l_save = n / 2;
    n_save = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( indx < 0 )
  {
    if ( indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = l_save;
      l_save = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k_save <= 1 )
    {
      if ( n_save == 1 )
      {
        i_save = 0;
        j_save = 0;
        indx = 0;
      }
      else
      {
        i_save = n_save;
        j_save = 1;
        n_save = n_save - 1;
        indx = 1;
      }
      i = i_save;
      j = j_save;
      return;
    }
    k_save = k_save - 1;
    l_save = k_save;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( indx == 1 )
  {
    l_save = k_save;
  }

  for ( ; ; )
  {

    i_save = 2 * l_save;

    if ( i_save == n_save )
    {
      j_save = l_save;
      l_save = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }
    else if ( i_save <= n_save )
    {
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k_save <= 1 )
    {
      break;
    }

    k_save = k_save - 1;
    l_save = k_save;
  }

  if ( n_save == 1 )
  {
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  }
  else
  {
    i_save = n_save;
    j_save = 1;
    n_save = n_save - 1;
    indx = 1;
    i = i_save;
    j = j_save;
  }

  return;
}
//****************************************************************************80

void sort_safe_rc ( int n, int &indx, int &i, int &j, int isgn, int &i_save,
  int &j_save, int &k_save, int &l_save, int &n_save )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_SAFE_RC externally ascending sorts a list of items.
//
//  Discussion:
//
//    This is a version of SORT_RC which does not rely on
//    storing certain work variables internally to the function.  This makes
//    the function somewhat more awkward to call, but easier to program
//    in a variety of languages, and safe to use in a parallel programming
//    environment, or in cases where the sorting of several vectors is to
//    be carried out at more or less the same time.
//
//    The actual list of data is not passed to the routine.  Hence this
//    routine may be used to sort integers, reals, numbers, names,
//    dates, shoe sizes, and so on.  After each call, the routine asks
//    the user to compare or interchange two items, until a special
//    return value signals that the sorting is completed.
//
//  Example:
//
//    n = 100;
//    indx = 0;
//
//    for ( ; ; )
//    {
//      sort_rc ( n, indx, i, j, isgn, i_save, j_save, k_save, 
//        l_save, n_save );
//
//      if ( indx < 0 )
//      {
//        isgn = 1;
//        if ( a[i-1] <= a[j-1] )
//        {
//          isgn = -1;
//        }
//      }
//      else if ( 0 < indx )
//      {
//        k      = a[i-1];
//        a[i-1] = a[j-1];
//        a[j-1] = k;
//      }
//      else
//      {
//        break;
//      }
//    }
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2015
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the length of the input list.
//
//    Input/output, int &INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int &I, &J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
//    Input/output, int &I_SAVE, &J_SAVE, &K_SAVE, &L_SAVE,
//    &N_SAVE, workspace needed by the routine.  Before calling the function,
//    the user should declare variables to hold these values, but should
//    not change them, and need not ever examine them.
//
{
//
//  INDX = 0: This is the first call.
//
  if ( indx == 0 )
  {
    i_save = 0;
    j_save = 0;
    k_save = n / 2;
    l_save = n / 2;
    n_save = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( indx < 0 )
  {
    if ( indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = l_save;
      l_save = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      indx = 2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k_save <= 1 )
    {
      if ( n_save == 1 )
      {
        i_save = 0;
        j_save = 0;
        indx = 0;
      }
      else
      {
        i_save = n_save;
        j_save = 1;
        n_save = n_save - 1;
        indx = 1;
      }
      i = i_save;
      j = j_save;
      return;
    }
    k_save = k_save - 1;
    l_save = k_save;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( indx == 1 )
  {
    l_save = k_save;
  }

  for ( ; ; )
  {
    i_save = 2 * l_save;

    if ( i_save == n_save )
    {
      j_save = l_save;
      l_save = i_save;
      indx = -1;
      i = i_save;
      j = j_save;
      return;
    }
    else if ( i_save <= n_save )
    {
      j_save = i_save + 1;
      indx = -2;
      i = i_save;
      j = j_save;
      return;
    }

    if ( k_save <= 1 )
    {
      break;
    }

    k_save = k_save - 1;
    l_save = k_save;
  }

  if ( n_save == 1 )
  {
    i_save = 0;
    j_save = 0;
    indx = 0;
    i = i_save;
    j = j_save;
  }
  else
  {
    i_save = n_save;
    j_save = 1;
    n_save = n_save - 1;
    indx = 1;
    i = i_save;
    j = j_save;
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
