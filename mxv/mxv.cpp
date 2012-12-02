# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
double cpu_time ( void );
double *matgen ( int m, int n );
double mxv_foriforj ( int m, int n, double a[], double x[], double y[] );
double mxv_forjfori ( int m, int n, double a[], double x[], double y[] );
void timestamp ( void );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MXV.
//
//  Discussion:
//
//    MXV computes a matrix-vector product in a number of ways, and reports
//    the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    mxv m n
//
//  Parameters:
//
//    Command line argument, int M, the number of rows in the matrix.
//
//    Command line argument, int N, the number of columns in the matrix.
//
{
  double *a;
  double cpu_seconds;
  int flop_count;
  int i;
  int m;
  double mflops;
  int n;
  double *x;
  double *y;

  timestamp ( );

  cout << "\n";
  cout << "MXV:\n";
  cout << "  C++ version\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
  cout << "\n";
  cout << "  Compute matrix vector products y = A*x.\n";
//
//  Get the number of rows, M.
//
  if ( 2 <= argc )
  {
    m = atoi ( argv[1] );
  } 
  else
  {
    cout << "\n";
    cout << "  Enter the number of rows, M\n";
    cin >> m;
  }
//
//  Get the number of columns, N.
//
  if ( 3 <= argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the number of columns, N\n";
    cin >> n;
  }
//
//  Record the amount of work.
//  Each of the M entries of Y requires N multiplies and N adds.
//
  flop_count = 2 * m * n;

  cout << "\n";
  cout << "  Number of matrix rows M =             " << m << "\n";
  cout << "  Number of matrix columns N =          " << n  << "\n";
  cout << "  Number of floating point operations = " << flop_count  << "\n";
//
//  Set A and X.
//
  a = matgen ( m, n );

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }

  y = new double[m];

  cout << "\n";
  cout << "  Method     Cpu Seconds       MegaFlopS\n";
  cout << "  ------  --------------  --------------\n";
//
//  FORIFORJ
//
  cpu_seconds = mxv_foriforj ( m, n, a, x, y );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = - 1.0;
  }

  cout << "  FORIFORJ"
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  FORJFORI
//
  cpu_seconds = mxv_forjfori ( m, n, a, x, y );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = - 1.0;
  }

  cout << "  FORJFORI"
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  Deallocate arrays.
//
  delete [] a;
  delete [] x;
  delete [] y;

  cout << "\n";
  cout << "MXV:\n";
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

double *matgen ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MATGEN generates a random matrix A and vector X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Output, double MATGEN[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;
  int k;
  int seed;

  a = new double[m*n];

  seed = 1325;
//
// Set the matrix A.
//
  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      seed = ( ( 3125 * seed ) % 65536 );
      a[k] = ( seed - 32768.0 ) / 16384.0;
      k = k + 1;
    }
  }
  return a;
}
//****************************************************************************80

double mxv_foriforj ( int m, int n, double a[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXV_FORIFORJ computes y = A * x, using FOR I, FOR J loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns
//    of the matrix.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double Y[M], the product vector.
//
//    Output, double MXV_FORIFORJ, the elapsed CPU time.
//
{
  double cpu_seconds;
  int i;
  int j;
  double time1;
  double time2;

  time1 = cpu_time ( );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  time2 = cpu_time ( );

  cpu_seconds = time2 - time1;

  return cpu_seconds;
}
//****************************************************************************80

double mxv_forjfori ( int m, int n, double a[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXV_FORJFORI computes y = A * x, using FOR J, FOR I loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns
//    of the matrix.
//
//    Input, double A[M*N], the matrix.
//
//    Input, double X[N], the vector to be multiplied.
//
//    Output, double Y[M], the product vector.
//
//    Output, double MXV_FORJFORI, the elapsed CPU time.
//
{
  double cpu_seconds;
  int i;
  int j;;
  double time1;
  double time2;

  time1 = cpu_time ( );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  time2 = cpu_time ( );

  cpu_seconds = time2 - time1;

  return cpu_seconds;
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
