# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include "wtime.hpp"

using namespace std;

int main ( );
void test01 ( );
int i4_power ( int i, int j );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WTIME_PRB.
//
//  Discussion:
//
//    WTIME_PRB tests the WTIME library.
//
//    CLOCK is a timing utility accessible to C codes. It returns the number of
//    "ticks" of processor time devoted to the user's job.  Dividing this by
//    CLOCKS_PER_SEC (the number of clock ticks per second) gives a rough value for
//    the elapsed CPU time.
//
//    In this program, we call a routine CPU_TIME, which conveniently converts
//    CLOCK's output to seconds for us.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "WTIME_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WTIME library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WTIME_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 times the RAND routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  double cpu1;
  double cpu2;
  int i;
  int i_rep;
  int n;
  int n_log;
  int n_log_min = 10;
  int n_log_max = 20;
  int n_max;
  int n_min;
  int n_rep = 5;
  double seconds;
  int *x;

  n_min = i4_power ( 2, n_log_min );
  n_max = i4_power ( 2, n_log_max );

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Time the RAND routine by computing N values.\n";
  cout << "  For a given N, repeat the computation 5 times.\n";
  cout << "\n";
  cout << "  Data vectors will be of minimum size " << n_min << "\n";
  cout << "  Data vectors will be of maximum size " << n_max << "\n";
  cout << "\n";
  cout << "  Times are measured in seconds.\n";
  cout << "\n";
  cout << "         N      Rep #1      Rep #2      Rep #2      Rep #4      Rep #5\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = i4_power ( 2, n_log );
    x = new int[n];

    cout << "  " << setw(8) << n;

    for ( i_rep = 0; i_rep < n_rep; i_rep++ )
    {
      seconds = wtime ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      seconds = wtime ( ) - seconds;

      cout << "  " << setw(10) << seconds;
    }
    cout << "\n";
    delete [] x;
  }

  return;
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
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
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
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
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
