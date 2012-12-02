# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <omp.h>

using namespace std;

int main ( );
void test01 ( int m, int n );
void matgen ( int m, int n, double a[], double x[] );
void mxv_plain ( int m, int n, double a[], double x[], double y[] );
void mxv_plain_openmp ( int m, int n, double a[], double x[], double y[] );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MXV_OPENMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  int n;

  timestamp ( );

  cout << "\n";
  cout << "MXV_OPENMP:\n";
  cout << "  C++/OpenMP version\n";
  cout << "\n";
  cout << "  Compute matrix vector products y = A*x.\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  cout << "\n";
  cout << "  Compare various algorithms:\n";
  cout << "\n";
  cout << "  MXV_PLAIN          - plain MxV coding.\n";
  cout << "  MXV_PLAIN_OPENMP  - plain MxV coding + OpenMP.\n";
  cout << "\n";
  cout << "  Algorithm                  M         N      Seconds\n";
//
//  N = M
//
  m = 10;

  for ( i = 1; i <= 3; i++ )
  {
    cout << "\n";

    n = m;
    test01 ( m, n );

    m = m * 10;
  }
//
//  N = 10 * M
//
  m = 1;

  for ( i = 1; i <= 4; i++ )
  {
    cout << "\n";

    n = 10 * m;
    test01 ( m, n );

    m = m * 10;
  }
//
//  M = 10 * N
//
  n = 1;

  for ( i = 1; i <= 4; i++ )
  {
    cout << "\n";

    m = 10 * n;
    test01 ( m, n );

    n = n * 10;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "MXV_OPENMP:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 compares various algorithms for a given matrix size MxN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2008
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
{
  double *a;
  double seconds;
  double *x;
  double *y;

  a = new double[m*n];
  x = new double[n];
  y = new double[m];

  matgen ( m, n, a, x );

  seconds = omp_get_wtime ( );
  mxv_plain ( m, n, a, x, y );
  seconds = omp_get_wtime ( ) - seconds;
  cout << "  MXV_PLAIN         "
    << "  " << setw(8) << m
    << "  " << setw(8) << n
    << "  " << setw(14) << seconds << "\n";

  seconds = omp_get_wtime ( );
  mxv_plain_openmp ( m, n, a, x, y );
  seconds = omp_get_wtime ( ) - seconds;
  cout << "  MXV_PLAIN_OPENMP "
    << "  " << setw(8) << m
    << "  " << setw(8) << n
    << "  " << setw(14) << seconds << "\n";

  delete [] a;
  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void matgen ( int m, int n, double a[], double x[] )

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
//    08 April 2008
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
//    Output, double A[M*N], the matrix.
//
//    Output, double X[N], the vector.
//
{
  int i;
  int j;
  int seed;

  seed = 1325;
//
// Set the matrix A.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      seed = ( 3125 * seed ) % 65536;
      a[i+j*m] = ( seed - 32768.0 ) / 16384.0;
    }
  }
//
//  Set X.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = i + 1;
  }

  return;
}
//****************************************************************************80

void mxv_plain ( int m, int n, double a[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXV_PLAIN computes y = A * x, using "plain" code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2008
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
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return;
}
//****************************************************************************80

void mxv_plain_openmp ( int m, int n, double a[], double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    MXV_PLAIN_OPENMP computes y = A * x, using OpenMP parallel directives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
//    Using OpenMP: Portable Shared Memory Parallel Processing,
//    MIT Press, 2007,
//    ISBN13: 978-0262533027,
//    LC: QA76.642.C49.
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
{
  int i;
  int j;

# pragma omp parallel shared ( m, n, a, x, y ) private ( i, j )
# pragma omp for
  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
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
