# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>

using namespace std;

# include "disk_grid.hpp"

//****************************************************************************80

double *disk_grid ( int n, double r, double c[2], int ng )

//****************************************************************************80
//
//  Purpose:
//
//    DISK_GRID computes grid points inside a disk.
//
//  Discussion:
//
//    The grid is defined by specifying the radius and center of the circle,
//    and the number of subintervals N into which the horizontal radius
//    should be divided.  Thus, a value of N = 2 will result in 5 points
//    along that horizontal line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R, the radius of the circle.
//
//    Input, double C[2], the coordinates of the center of the circle.
//
//    Input, int NG, the number of grid points, as determined by
//    DISK_GRID_COUNT.
//
//    Output, double DISK_GRID[2*NG], the grid points inside the circle.
//
{
  double *cg;
  int i;
  int j;
  int p;
  double x;
  double y;

  cg = new double[2*ng];

  p = 0;

  for ( j = 0; j <= n; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );

    cg[0+2*p] = x;
    cg[1+2*p] = y;
    p = p + 1;

    if ( 0 < j )
    {
      cg[0+2*p] = x;
      cg[1+2*p] = 2.0 * c[1] - y;
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );

      if ( r * r < pow ( x - c[0], 2 ) + pow ( y - c[1], 2 ) )
      {
        break;
      }

      cg[0+2*p] = x;
      cg[1+2*p] = y;
      p = p + 1;
      cg[0+2*p] = 2.0 * c[0] - x;
      cg[1+2*p] = y;
      p = p + 1;

      if ( 0 < j )
      {
        cg[0+2*p] = x;
        cg[1+2*p] = 2.0 * c[1] - y;
        p = p + 1;
        cg[0+2*p] = 2.0 * c[0] - x;
        cg[1+2*p] = 2.0 * c[1] - y;
        p = p + 1;
      }
    }
  }
  return cg;
}
//****************************************************************************80

int disk_grid_count ( int n, double r, double c[2] )

//****************************************************************************80
//
//  Purpose:
//
//    DISK_GRID_COUNT counts the grid points inside a disk.
//
//  Discussion:
//
//    The grid is defined by specifying the radius and center of the circle,
//    and the number of subintervals N into which the horizontal radius
//    should be divided.  Thus, a value of N = 2 will result in 5 points
//    along that horizontal line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R, the radius of the circle.
//
//    Input, double C[2], the coordinates of the center of the circle.
//
//    Output, int DISK_GRID_COUNT, the number of grid points inside 
//    the circle.
//
{
  int i;
  int j;
  int ng;
  double x;
  double y;

  ng = 0;

  for ( j = 0; j <= n; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );
    ng = ng + 1;

    if ( 0 < j )
    {
      ng = ng + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );

      if ( r * r < pow ( x - c[0], 2 ) + pow ( y - c[1], 2 ) )
      {
        break;
      }

      ng = ng + 1;
      ng = ng + 1;
      if ( 0 < j )
      {
        ng = ng + 1;
        ng = ng + 1;
      }
    }
  }
  return ng;
}
//****************************************************************************80

double *disk_grid_fibonacci ( int n, double r, double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    DISK_GRID_FIBONACCI computes Fibonacci grid points inside a disk.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Swinbank, James Purser,
//    Fibonacci grids: A novel approach to global modelling,
//    Quarterly Journal of the Royal Meteorological Society,
//    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
//
//  Parameters:
//
//    Input, int N, the number of points desired.
//
//    Input, double R, the radius of the circle.
//
//    Input, double C[2], the coordinates of the center of the circle.
//
//    Output, double DISK_GRID_FIBONACCI[2*N], the grid points.
//
{
  double *g;
  double gr;
  double gt;
  int i;
  double phi;
  const double pi = 3.141592653589793;
  double r0;

  r0 = r / sqrt ( ( double ) ( n ) - 0.5 );
  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  g = new double[2*n];

  for ( i = 0; i < n; i++ )
  {
    gr = r0 * sqrt ( ( double ) ( i + 1 ) - 0.5 );
    gt = 2.0 * pi * ( double ) ( i + 1 ) / phi;
    g[0+i*2] = c[0] + gr * cos ( gt );
    g[1+i*2] = c[1] + gr * sin ( gt );
  }

  return g;
}
//****************************************************************************80

void r82vec_print_part ( int n, double a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PRINT_PART prints "part" of an R82VEC.
//
//  Discussion:
//
//    The user specifies MAX_PRINT, the maximum number of lines to print.
//
//    If N, the size of the vector, is no more than MAX_PRINT, then
//    the entire vector is printed, one entry per line.
//
//    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
//    followed by a line of periods suggesting an omission,
//    and the last entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[2*N], the vector to be printed.
//
//    Input, int MAX_PRINT, the maximum number of lines
//    to print.
//
//    Input, string TITLE, a title.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(8) << i
           << "  " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2]  << "\n";
    }
    cout << "  ........  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*2]
         << "  " << setw(14) << a[1+i*2]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*2]
           << "  " << setw(14) << a[1+i*2]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*2]
         << "  " << setw(14) << a[1+i*2] 
         << "  " << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
