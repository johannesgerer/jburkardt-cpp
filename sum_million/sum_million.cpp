# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( void );
double cpu_time ( void );
double *set_up ( int n );
void sum_up ( int n, double x[], double *total, double *ctime );
void timestamp ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUM_MILLION.
//
//  Discussion:
//
//    This code estimates the power of a computer by summing the integers
//    from 1 to 1,000,000.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  double error;
  double exact = 500000500000.0;
  int i;
  double mflops;
  int n = 1000000;
  double total;
  double *x;

  cout << "\n";
  timestamp ( );
  cout << "\n";
  cout << "SUM_MILLION\n";
  cout << "  C++ version\n";
  cout << "  Sum the integers from 1 to 1,000,000.\n";
  cout << "  Correct answer is 500000500000.\n";

  x = set_up ( n );

  cout << "\n";
  cout << "         N      CPU time        MFLOPS          ERROR\n";
  cout << "                (seconds)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    sum_up ( n, x, &total, &ctime );

    mflops = ( double ) ( n ) / 1000000.0 / ctime;

    error = total - exact;

    cout << "  " << setw(8) << n
         << "  " << setw(14) << ctime
         << "  " << setw(14) << mflops 
         << "  " << setw(14) << error << "\n";
  }

  delete [] x;

  cout << "\n";
  cout << "SUM_MILLION:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double cpu_time ( void )

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

double *set_up ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    SET_UP sets up the data for the SUM_MILLION program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to define.
//
//    Output, double X[N], a vector which contains the values 1 through N.
//
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  return x;
}
//****************************************************************************80

void sum_up ( int n, double x[], double *total, double *ctime )

//****************************************************************************80
//
//  Purpose:
//
//    SUM_UP carries out the sum for the SUM_MILLION program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values to define.
//
//    Input, double X[N], the data to be summed.
//
//    Output, double TOTAL, the sum of the data.
//
//    Output, double CTIME, the cpu time required to sum the data.
//
{
  double ctime1;
  double ctime2;
  int i;

  ctime1 = cpu_time ( );

  *total = 0.0;
  for ( i = 0; i < n; i++ )
  {
    *total = *total + x[i];
  }

  ctime2 = cpu_time ( );

  *ctime = ctime2 - ctime1;

  return;
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
