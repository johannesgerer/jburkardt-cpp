# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <omp.h>

using namespace std;

int main ( void );
void test01 ( void );
int i4_power ( int i, int j );
void timestamp ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TIMER_OMP_GET_WTIME.
//
//  Discussion:
//
//    TIMER_OMP_GET_WTIME uses OMP_GET_WTIME as the timer.
//
//    OMP_GET_WTIME is a timing utility accessible to C codes that
//    support OpenMP.  It returns the elapsed wallclock time in seconds.
//
//    Here, we run on as many threads as there are processors.  We could
//    force the number of threads to be 1 to make a better comparison to
//    timers that run on a single processor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  int proc_num;
  int thread_num;

  timestamp ( );

  cout << "\n";
  cout << "TIMER_OMP_GET_WTIME\n";
  cout << "  C++ version\n";
  cout << "\n"; 
  cout << "  Demonstrate the use of the OMP_GET_WTIME timer.\n";
  cout << "\n";
  cout << "  omp_get_wtime ( ) is an OpenMP library function.\n";
  cout << "\n";
  cout << "  It returns the elapsed wall clock time in seconds.\n";
//
//  How many processors are available?
//
  proc_num = omp_get_num_procs ( );

  cout << "\n";
  cout << "  The number of processors available:\n";
  cout << "  OMP_GET_NUM_PROCS () = " << proc_num << "\n";
  
  thread_num = proc_num;

  cout << "\n";
  cout << "  OMP_SET_NUM_THREADS requests " << thread_num << " threads.\n";

  omp_set_num_threads ( thread_num );

  test01 ( );

  cout << "\n";
  cout << "TIMER_CLOCK\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 times the C RAND routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int i_rep;
  int n;
  int n_log;
  int n_log_min = 10;
  int n_log_max = 20;
  int n_max;
  int n_min;
  int n_rep = 5;
  double wtime1;
  double wtime2;
  int *x;

  n_min = i4_power ( 2, n_log_min );
  n_max = i4_power ( 2, n_log_max );

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Time the C RAND routine by computing N values.\n";
  cout << "  For a given N, repeat the computation 5 times.\n";
  cout << "\n";
  cout << "  Data vectors will be of minimum size " << n_min << "\n";
  cout << "  Data vectors will be of maximum size " << n_max << "\n";
  cout << "\n";
  cout << "  Wall clock times are in seconds.\n";
  cout << "\n";
  cout << "         N       Rep #1       Rep #2       Rep #2       Rep #4       Rep #5\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    n = i4_power ( 2, n_log );
    x = new int[n];

    cout << "  " << setw(8) << n;

    for ( i_rep = 0; i_rep < n_rep; i_rep++ )
    {
      wtime1 = omp_get_wtime ( );

      for ( i = 0; i < n; i++ )
      {
        x[i] = rand ( );
      }

      wtime2 = omp_get_wtime ( );

      cout << "  " << setw(11) << wtime2 - wtime1;
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

void timestamp ( void )

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
