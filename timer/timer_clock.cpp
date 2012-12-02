# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( );
void test01 ( );
double cpu_time ( );
int i4_power ( int i, int j );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TIMER_CLOCK.
//
//  Discussion:
//
//    TIMER_CLOCK uses CLOCK as the timer.
//
//    CLOCK is a timing utility accessible to C++ codes. It returns the number of
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
  cout << "TIMER_CLOCK\n";
  cout << "  C++ version\n";
  cout << "  Demonstrate the use of the CLOCK timer.\n";
  cout << "\n";
  cout << "  CLOCK is a  C++ library routine\n";
  cout << "  (defined in <ctime>)\n";
  cout << "\n";
  cout << "  It returns the processor time used by the program\n";
  cout << "  since the beginning of program execution.\n";
  cout << "  Divide this by CLOCKS_PER_SEC to convert to seconds.\n";
  cout << "\n";
  cout << "  CLOCK is a crude timer, and results less than\n";
  cout << "  a tenth of a second are probably not reliable.\n";
  cout << "\n";
  cout << "  The number of clock ticks per second is " 
       << CLOCKS_PER_SEC << "\n";

  test01 ( );

  cout << "\n";
  cout << "TIMER_CLOCK\n";
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
//    16 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cpu;
  double cpu_ave;
  double cpu_max;
  double cpu_min;
  int i;
  int n;
  int n_log;
  int n_log_min = 0;
  int n_log_max = 20;
  int n_max;
  int n_min;
  int rep;
  int rep_num = 5;
  int *x;

  n_min = i4_power ( 2, n_log_min );
  n_max = i4_power ( 2, n_log_max );

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Time the RAND function by computing N values.\n";
  cout << "  For a given N, repeat the computation 5 times.\n";
  cout << "\n";
  cout << "  Data vectors will be of minimum size " << n_min << "\n";
  cout << "  Data vectors will be of maximum size " << n_max << "\n";
  cout << "\n";
  cout << "  CPU times are in seconds.\n";
  cout << "\n";
  cout << "         N      Rep #1      Rep #2      Rep #2      Rep #4      Rep #5"
       << "         Min         Ave         Max\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = i4_power ( 2, n_log );
    x = new int[n];

    cout << "  " << setw(8) << n;

    cpu_min = r8_huge ( );
    cpu_max = 0.0;
    cpu_ave = 0.0;

    for ( rep = 1; rep <= rep_num; rep++ )
    {
      cpu = cpu_time ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      cpu = cpu_time ( ) - cpu;
      cpu_max = r8_max ( cpu_max, cpu );
      cpu_min = r8_min ( cpu_min, cpu );
      cpu_ave = cpu_ave + cpu;

      cout << "  " << setw(10) << cpu;
    }
    cpu_ave = cpu_ave / ( double ) rep_num;

    cout << "  " << setw(10) << cpu_min
         << "  " << setw(10) << cpu_ave
         << "  " << setw(10) << cpu_max << "\n";

    delete [] x;
  }

  return;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Discussion:
//
//    The data available to this routine through "CLOCK" is not very reliable,
//    and hence the values of CPU_TIME returned should not be taken too 
//    seriously, especially when short intervals are being timed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) 
        / ( double ) CLOCKS_PER_SEC;

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

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
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
