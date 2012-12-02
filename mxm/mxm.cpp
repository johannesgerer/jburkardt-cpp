# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
double cpu_time ( );
double *matgen ( int m, int n, int *seed );
double mxm_ijk ( int n1, int n2, int n3, double b[], double c[] );
double mxm_ikj ( int n1, int n2, int n3, double b[], double c[] );
double mxm_jik ( int n1, int n2, int n3, double b[], double c[] );
double mxm_jki ( int n1, int n2, int n3, double b[], double c[] );
double mxm_kij ( int n1, int n2, int n3, double b[], double c[] );
double mxm_kji ( int n1, int n2, int n3, double b[], double c[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MXM.
//
//  Discussion:
//
//    MXV computes a matrix-matrix product in a number of ways, and reports
//    the elapsed CPU time.
//
//    The multiplication carried out is
//
//      A(1:N1,1:N3) = B(1:N1,1:N2) * C(1:N2,1:N3)
//
//    where B and C are real double precision matrices whose entries
//    are assigned randomly.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    mxm n1 n2 n3
//
//  Parameters:
//
//    Command line argument, int N1, N2, N3, defines the number of
//    rows and columns in the two matrices.
//
{
  double *b;
  double *c;
  double cpu_seconds;
  int flop_count;
  double mflops;
  int n1;
  int n2;
  int n3;
  int seed;
  double time_estimate;

  timestamp ( );

  cout << "\n";
  cout << "MXM:\n";
  cout << "  C++ version\n";
  cout << "  Compute matrix-matrix product A = B * C\n";
//
//  Get N1.
//
  if ( 1 < argc )
  {
    n1 = atoi ( argv[1] );
  } 
  else
  {
    cout << "\n";
    cout << "  Enter N1, the number of rows in B.\n";
    cin >> n1;
  }
//
//  Get N2.
//
  if ( 2 < argc )
  {
    n2 = atoi ( argv[2] );
  } 
  else
  {
    cout << "\n";
    cout << "  Enter N2, the number of columns in B and rows in C.\n";
    cin >> n2;
  }
//
//  Get N3.
//
  if ( 3 < argc )
  {
    n3 = atoi ( argv[3] );
  } 
  else
  {
    cout << "\n";
    cout << "  Enter N3, the number of columns in C.\n";
    cin >> n3;
  }
//
//  Record the amount of work.
//  Each of the N1 * N3 entries of A requires N2 multiplies and N2 adds.
//
  flop_count = 2 * n1 * n2 * n3;

  cout << "\n";
  cout << "  Matrix B is " << n1 << " by " << n2 << "\n";
  cout << "  Matrix C is " << n2 << " by " << n3 << "\n";
  cout << "  Matrix A will be " << n1 << " by " << n3 << "\n";
  cout << "\n";
  cout << "  Number of floating point operations = " << flop_count << "\n";
  time_estimate = ( double ) ( flop_count ) / 2.6E+09;
  cout << "  Estimated CPU time is " << time_estimate << " seconds.\n";
//
//  Set B and C.
//
  seed = 1325;
  b = matgen ( n1, n2, &seed );
  c = matgen ( n2, n3, &seed );
  cout << "\n";
  cout << "  Method     Cpu Seconds       MegaFlopS\n";
  cout << "  ------  --------------  --------------\n";
//
//  IJK
//
  cpu_seconds = mxm_ijk ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  IJK   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  IKJ
//
  cpu_seconds = mxm_ikj ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  IKJ   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  JIK
//
  cpu_seconds = mxm_jik ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  JIK   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  JKI
//
  cpu_seconds = mxm_jki ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  JKI   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  KIJ
//
  cpu_seconds = mxm_kij ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  KIJ   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  KJI
//
  cpu_seconds = mxm_kji ( n1, n2, n3, b, c );

  if ( 0.0 < cpu_seconds )
  {
    mflops = ( double ) ( flop_count ) / cpu_seconds / 1000000.0;
  }
  else
  {
    mflops = -1.0;
  }

  cout << "  KJI   "
       << "  " << setw(14) << cpu_seconds
       << "  " << setw(14) << mflops << "\n";
//
//  Deallocate arrays.
//
  delete [] b;
  delete [] c;
//
//  Terminate.
//
  cout << "\n";
  cout << "MXM:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
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

double *matgen ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MATGEN generates a random matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2010
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
//    Input, int *SEED, a seed for the random number
//    generator.
//
//    Output, double MATGEN[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];
//
//  Set the matrix A.
//
  for ( j = 0;j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      *seed = ( ( 3125 * *seed ) % 65536 );
      a[i+j*m] = ( *seed - 32768.0 ) / 16384.0;
    }
  }

  return a;
}
//****************************************************************************80

double mxm_ijk ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_IJK computes A = B * C using FOR I, FOR J, FOR K loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_IJK, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      for ( k = 0; k < n2; k++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
}
//****************************************************************************80

double mxm_ikj ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_IKJ computes A = B * C using FOR I, FOR K, FOR J loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_IKJ, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( i = 0; i < n1; i++ )
  {
    for ( k = 0; k < n2; k++ )
    {
      for ( j = 0; j < n3; j++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
}
//****************************************************************************80

double mxm_jik ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_JIK computes A = B * C using FOR J, FOR I, FOR K loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_JIK, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      for ( k = 0; k < n2; k++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
}
//****************************************************************************80

double mxm_jki ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_JKI computes A = B * C using FOR J, FOR K, FOR I loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_JKI, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( j = 0; j < n3; j++ )
  {
    for ( k = 0; k < n2; k++ )
    {
      for ( i = 0; i < n1; i++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
}
//****************************************************************************80

double mxm_kij ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_KIJ computes A = B * C using FOR K, FOR I, FOR J loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_KIJ, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( k = 0; k < n2; k++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      for ( j = 0; j < n3; j++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
}
//****************************************************************************80

double mxm_kji ( int n1, int n2, int n3, double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXM_KJI computes A = B * C using FOR K, FOR J, FOR I loops.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, define the orders of the
//    matrices.
//
//    Input, double B[N1*N2], C[N2*N3], the factor matrices.
//
//    Output, double MXM_KJI, the elapsed CPU time.
//
{
  double *a;
  double cpu_seconds;
  int i;
  int j;
  int k;

  a = new double[n1*n3];

  for ( j = 0; j < n3; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      a[i+j*n1] = 0.0;
    }
  }

  cpu_seconds = cpu_time ( );

  for ( k = 0; k < n2; k++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      for ( i = 0; i < n1; i++ )
      {
        a[i+j*n1] = a[i+j*n1] + b[i+k*n1] * c[k+j*n2];
      }
    }
  }

  cpu_seconds = cpu_time ( ) - cpu_seconds;

  delete [] a;

  return cpu_seconds;
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
