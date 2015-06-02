# include <cstdlib>
# include <iostream>
# include <ctime>

using namespace std;

int main ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CELLULAR_AUTOMATON.
//
//  Discussion:
//
//    This program carries out iterations of the 1D cellular automaton
//    known as rule 30.
//
//    Given an initial linear array of 0's and 1's, rule 30 produces a new
//    array using the rules:
//
//      111  110  101  100  011  010  001  000
//       V    V    V    V    V    V    V    V
//       0    0    0    1    1    1    1    0     
//
//    Note that there are 256 = 2^8 possible ways to fill in this output
//    chart, and that rule 30 gets its index by the fact that
//    (0,0,0,1,1,1,1,0) can be interpreted as the binary representation of 30.
//
//    For instance, if the current values of X(4), X(5) and X(6) are
//    0, 1 and 1, respectively, then the new value of X(5) will be 1.
//
//    The first and last entries of the array must be treated specially, since
//    they don't have a left or right neighbor.  One simple treatment is 
//    to assume that there are phantom neighbors whose values are both 0.
//    Another is to enforce periodic boundary conditions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    A New Kind of Science,
//    Wolfram Media, 2002,
//    ISBN13: 978-1579550080,
//    LC: QA267.5.C45.W67.
//
{
  int i;
  int j;
  int n;
  int step_num;
  char *x;
  char *x_old;

  timestamp ( );
  cout << "\n";
  cout << "CELLULAR_AUTOMATON:\n";
  cout << "  C++ version.\n";

  n = 80;
  step_num = 80;

  x = new char[n+2];
  x_old = new char[n+2];

  for ( i = 0; i <= n + 1; i++ )
  {
    x[i] = ' ';
  }
  x[40] = '*';

  for ( i = 1; i <= n; i++ )
  {
    cout << x[i];
  }
  cout << "\n";

  for ( j = 1; j <= step_num; j++ )
  {
    for ( i = 0; i < n + 2; i++ )
    {
      x_old[i] = x[i];
    }
    for ( i = 1; i <= n; i++ )
    {
//
//  The transformation rules are:
//
//  111  110  101  100  011  010  001  000
//   |    |    |    |    |    |    |    |
//   0    0    0    1    1    1    1    0
//
//  which means this rule has binary code 00011110 = 16 + 8 + 4 + 2 = 30
//
      if ( ( x_old[i-1] == ' ' && x_old[i] == ' ' && x_old[i+1] == '*' ) ||
           ( x_old[i-1] == ' ' && x_old[i] == '*' && x_old[i+1] == ' ' ) ||
           ( x_old[i-1] == ' ' && x_old[i] == '*' && x_old[i+1] == '*' ) ||
           ( x_old[i-1] == '*' && x_old[i] == ' ' && x_old[i+1] == ' ' ) )
      {
        x[i] = '*';
      }
      else
      {
        x[i] = ' ';
      }
    }
//
//  Enforce periodic boundary conditions.
//
    x[0] = x[n];
    x[n+1] = x[1];

    for ( i = 1; i <= n; i++ )
    {
      cout << x[i];
    }
    cout << "\n";
  }
//
//  Free memory.
//
  delete [] x;
  delete [] x_old;
//
//  Terminate.
//
  cout << "\n";
  cout << "CELLULAR_AUTOMATON:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
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
