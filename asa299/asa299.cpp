# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "asa299.hpp"

//****************************************************************************80

void simplex_lattice_point_next ( int n, int t, bool *more, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SIMPLEX_LATTICE_POINT_NEXT generates lattice points in a simplex.
//
//  Discussion:
//
//    The simplex is defined by N-dimensional points X such that:
//
//        0 <= X(1:N)
//
//    and
//
//      sum ( X(1:N) ) <= T
//
//    where T is an integer.
//
//    Lattice points are points X which satisfy the simplex conditions and
//    for which all the components are integers.
//
//    This routine generates all the lattice points in a given simplex, one at 
//    a time, in a reverse lexicographic order.
//
//    To use the routine, initialize by setting N and T to appropriate values, 
//    and MORE to FALSE.  The initial value of X is not important.
//
//    Call the routine. On return, X will contain the first lattice point in 
//    the simplex.  If MORE is TRUE, then the routine may be called again to 
//    get the next point.  In fact, as long as the output value of MORE is 
//    TRUE, there is at least one more lattice point that can be found by 
//    making another call.  When MORE is returned as FALSE, then there are no 
//    more lattice points; the value of X returned at that time is the 
//    "last" such point.
//
//    During the computation of a sequence of lattice points, the user should 
//    not change the values of N, T, MORE or X.  
//
//    The output for N = 3, T = 4 would be:
//
//       1    4  0  0
//       2    3  1  0
//       3    3  0  1
//       4    2  2  0
//       5    2  1  1
//       6    2  0  2
//       7    1  3  0
//       8    1  2  1
//       9    1  1  2
//      10    1  0  3
//      11    0  4  0
//      12    0  3  1
//      13    0  2  2
//      14    0  1  3
//      15    0  0  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Scott Chasalow, Richard Brand,
//    Algorithm AS 299:
//    Generation of Simplex Lattice Points,
//    Applied Statistics,
//    Volume 44, Number 4, 1995, pages 534-545.
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
//    Input, int N, the spatial dimension.
//    N must be positive.
//
//    Input, int T, the characteristic of the simplex.
//    T must be nonnegative.
//
//    Input/output, int *MORE, initialized to FALSE by the user to
//    begin a sequence of calculations, returned by the routine as TRUE,
//    if there are more values of X that can be calculated, or FALSE
//    if the accompanying value of X is the last one for this sequence.
//
//    Input/output, int X[N], not initialized by the user, but not
//    changed by the user on subsequent calls.  The routine returns
//    a new point on each call, and on subsequent calls uses the input
//    value (old point) to compute the output value (next point).
//
{
  int i;
  int j;

  if ( !(*more) )
  {
    if ( n < 1 )
    {
      cout << "\n";
      cout << "SIMPLEX_LATTICE_POINT_NEXT - Fatal error!\n";
      cout << "  N < 1.\n";
      exit ( 1 );
    }

    if ( t < 0 )
    {
      cout << "\n";
      cout << "SIMPLEX_LATTICE_POINT_NEXT - Fatal error!\n";
      cout << "  T < 0.\n";
      exit ( 1 );
    }

    *more = 1;
    j = 1;

    x[0] = t;
    for ( i = 1; i < n; i++ )
    {
      x[i] = 0;
    }
//
//  The first point can actually also be the last!
//
    if ( n == 1 )
    {
      *more = 0;
    }
  }
  else
  {
//
//  Search X(N-1 down to 1) for the first nonzero element.
//  If none, then terminate.  (This should not happen!)
//  Otherwise, set J to this index.
//  Decrement X(J) by 1.
//  Set X(J+1:N) to (T-X(1:J),0,0,...0).
//
    j = n - 1;

    for ( i = n - 2; 0 <= i; i-- )
    {
      if ( 0 < x[i] )
      {
        j = i;
        break;
      }
    }

    if ( j == n - 1 )
    {
      cout << "\n";
      cout << "SIMPLEX_LATTICE_POINT_NEXT - Fatal error!\n";
      cout << "  The input X vector is nonpositive in all entries\n";
      cout << "  except possibly the last one.\n";
      cout << "\n";
      cout << "  Perhaps the user has miscalled the routine\n";
      cout << "  or altered data between calls.\n";
      cout << "\n";
      cout << "ABNORMAL TERMINATION.\n";
      exit ( 1 );
    }

    x[j] = x[j] - 1;
    x[j+1] = t;
    for ( i = 0; i <= j; i++ )
    {
      x[j+1] = x[j+1] - x[i];
    }
    for ( i = j+2; i < n; i++ )
    {
      x[i] = 0;
    } 
//
//  Is this the last point?
//
    if ( x[n-1] == t )
    {
      *more = 0;
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
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
