# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

# include <omp.h>

using namespace std;

int main ( );
int i4_min ( int i1, int i2 );
void i4pp_delete ( int **a, int m, int n );
int **i4pp_new ( int m, int n );
void timestamp ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose
//
//    MAIN is the main program for MANDELBROT_OPENMP.
//
//  Discussion:
//
//    MANDELBROT_OPENMP computes an image of the Mandelbrot set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    Local, int COUNT_MAX, the maximum number of iterations taken
//    for a particular pixel.
//
{
  int m = 500;
  int n = 500;

  int **b;
  int c;
  int c_max;
  int **count;
  int count_max = 2000;
  int **g;
  int i;
  int ierror;
  int j;
  int jhi;
  int jlo;
  int k;
  string filename = "mandelbrot.ppm";
  ofstream output;
  int **r;
  double wtime;
  double wtime_total;
  double x_max =   1.25;
  double x_min = - 2.25;
  double x;
  double x1;
  double x2;
  double y_max =   1.75;
  double y_min = - 1.75;
  double y;
  double y1;
  double y2;

  b = i4pp_new ( m, n );
  count = i4pp_new ( m, n );
  g = i4pp_new ( m, n );
  r = i4pp_new ( m, n );

  timestamp ( );
  cout << "\n";
  cout << "MANDELBROT_OPENMP\n";
  cout << "  C++/OpenMP version\n";
  cout << "\n";
  cout << "  Create an ASCII PPM image of the Mandelbrot set.\n";
  cout << "\n";
  cout << "  For each point C = X + i*Y\n";
  cout << "  with X range [" << x_min << "," << x_max << "]\n";
  cout << "  and  Y range [" << y_min << "," << y_max << "]\n";
  cout << "  carry out " << count_max << " iterations of the map\n";
  cout << "  Z(n+1) = Z(n)^2 + C.\n";
  cout << "  If the iterates stay bounded (norm less than 2)\n";
  cout << "  then C is taken to be a member of the set.\n";
  cout << "\n";
  cout << "  An ASCII PPM image of the set is created using\n";
  cout << "    M = " << m << " pixels in the X direction and\n";
  cout << "    N = " << n << " pixels in the Y direction.\n";

  wtime = omp_get_wtime ( );
/*
  Carry out the iteration for each pixel, determining COUNT.
*/
# pragma omp parallel \
  shared ( b, count, count_max, g, r, x_max, x_min, y_max, y_min ) \
  private ( i, j, k, x, x1, x2, y, y1, y2 )
{
# pragma omp for

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x = ( ( double ) (     j - 1 ) * x_max   
          + ( double ) ( m - j     ) * x_min ) 
          / ( double ) ( m     - 1 );

      y = ( ( double ) (     i - 1 ) * y_max   
          + ( double ) ( n - i     ) * y_min ) 
          / ( double ) ( n     - 1 );

      count[i][j] = 0;

      x1 = x;
      y1 = y;

      for ( k = 1; k <= count_max; k++ )
      {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2 * x1 * y1 + y;

        if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
        {
          count[i][j] = k;
          break;
        }
        x1 = x2;
        y1 = y2;
      }

      if ( ( count[i][j] % 2 ) == 1 )
      {
        r[i][j] = 255;
        g[i][j] = 255;
        b[i][j] = 255;
      }
      else
      {
        c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt ( 
          ( ( double ) ( count[i][j] ) / ( double ) ( count_max ) ) ) ) ) );
        r[i][j] = 3 * c / 5;
        g[i][j] = 3 * c / 5;
        b[i][j] = c;
      }
    }
  }
}

  wtime = omp_get_wtime ( ) - wtime;
  cout << "\n";
  cout << "  Time = " << wtime << " seconds.\n";
/*
  Write data to an ASCII PPM file.
*/
  output.open ( filename.c_str ( ) );

  output << "P3\n";
  output << n << "  " << m << "\n";
  output << 255 << "\n";
  for ( i = 0; i < m; i++ )
  {
    for ( jlo = 0; jlo < n; jlo = jlo + 4 )
    {
      jhi = i4_min ( jlo + 4, n );
      for ( j = jlo; j < jhi; j++ )
      {
        output << "  " << r[i][j]
               << "  " << g[i][j]
               << "  " << b[i][j] << "\n";
      }
      output << "\n";
    }
  }

  output.close ( );
  cout << "\n";
  cout << "  Graphics data written to \"" << filename << "\".\n";
/*
  Free memory.
*/
  i4pp_delete ( b, m, n );
  i4pp_delete ( count, m, n );
  i4pp_delete ( g, m, n );
  i4pp_delete ( r, m, n );
/*
  Terminate.
*/
  cout << "\n";
  cout << "MANDELBROT_OPENMP\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void i4pp_delete ( int **a, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4PP_DELETE frees memory associated with an I4PP.
//
//  Discussion:
//
//    This function releases the memory associated with an array that was 
//    created by a command like
//      int **a;
//      a = i4pp_new ( m, n );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int **A, a pointer to the pointers to the array.
//
//    Input, int M, N, the number of rows and columns in the array.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    delete [] a[i];
  }

  delete [] a;

  return;
}
//****************************************************************************80

int **i4pp_new ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4PP_NEW allocates a new I4PP.
//
//  Discussion:
//
//    A declaration of the form
//      int **a;
//    is necesary.  Then an assignment of the form:
//      a = i4pp_new ( m, n );
//    allows the user to assign entries to the matrix using typical
//    2D array notation:
//      a[2][3] = 17;
//      y = a[1][0];
//    and so on.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the matrix.
//
//    Output, int **I4PP_NEW, a pointer to the pointers to the array.
//
{
  int **a;
  int i;

  a = new int *[m];

  if ( a == NULL )
  {
    cerr << "\n";
    cerr << "I4PP_NEW - Fatal error!\n";
    cerr << "  Unable to allocate row pointer array.\n";
    exit ( 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    a[i] = new int[n];
    if ( a[i] == NULL )
    {
      cerr << "\n";
      cerr << "I4PP_NEW - Fatal error!\n";
      cerr << "  Unable to allocate row array.\n";
      exit ( 1 );
    }
  }

  return a;
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
