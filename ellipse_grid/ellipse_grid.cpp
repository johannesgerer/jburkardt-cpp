# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <fstream>

using namespace std;

# include "ellipse_grid.hpp"

//****************************************************************************80

double *ellipse_grid ( int n, double r[2], double c[2], int ng )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID generates grid points inside an ellipse.
//
//  Discussion:
//
//    The ellipse is specified as
//
//      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
//
//    The user supplies a number N.  There will be N+1 grid points along
//    the shorter axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R[2], the half axis lengths.
//
//    Input, double C[2], the center of the ellipse.
//
//    Input, int NG, the number of grid points inside the ellipse.
//
//    Output, double ELLIPSE_GRID[2*NG], the grid points.
//
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double *xy;
  double y;

  xy = new double[2*ng];

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = i4_ceiling ( r[0] / r[1] ) * ( double ) ( n );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    xy[0+p*2] = x;
    xy[1+p*2] = y;
    p = p + 1;

    if ( 0 < j )
    {
      xy[0+p*2] = x;
      xy[1+p*2] = 2.0 * c[1] - y;
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      xy[0+p*2] = x;
      xy[1+p*2] = y;
      p = p + 1;
      xy[0+p*2] = 2.0 * c[0] - x;
      xy[1+p*2] = y;
      p = p + 1;

      if ( 0 < j )
      {
        xy[0+p*2] = x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
        xy[0+p*2] = 2.0 * c[0] - x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
      }
    }
  }
  return xy;
}
//****************************************************************************80

int ellipse_grid_count ( int n, double r[2], double c[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
//
//  Discussion:
//
//    The ellipse is specified as
//
//      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
//
//    The user supplies a number N.  There will be N+1 grid points along
//    the shorter axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R[2], the half axis lengths.
//
//    Input, double C[2], the center of the ellipse.
//
//    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside 
//    the ellipse.
//
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double y;

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = i4_ceiling ( r[1] / r[0] ) * ( double ) ( n );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = i4_ceiling ( r[0] / r[1] ) * ( double ) ( n );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    p = p + 1;

    if ( 0 < j )
    {
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 ) 
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      p = p + 1;
      p = p + 1;

      if ( 0 < j )
      {
        p = p + 1;
        p = p + 1;
      }
    }
  }

  return p;
}
//****************************************************************************80

int i4_ceiling ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CEILING rounds an R8 up to the next I4.
//
//  Example:
//
//    X        I4_CEILING(X)
//
//   -1.1      -1
//   -1.0      -1
//   -0.9       0
//   -0.1       0
//    0.0       0
//    0.1       1
//    0.9       1
//    1.0       1
//    1.1       2
//    2.9       3
//    3.0       3
//    3.14159   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose ceiling is desired.
//
//    Output, int I4_CEILING, the ceiling of X.
//
{
  int value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1;
  }

  return value;
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
