# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

# include "ball_grid.hpp"

//****************************************************************************80

double *ball_grid ( int n, double r, double c[3], int ng )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_GRID computes grid points inside a ball.
//
//  Discussion:
//
//    The grid is defined by specifying the radius and center of the ball,
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
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R, the radius of the ball.
//
//    Input, double C[3], the coordinates of the center of the ball.
//
//    Input, int NG, the number of grid points, as determined by
//    BALL_GRID_COUNT.
//
//    Output, double BALL_GRID[3*NG], the grid points inside the ball.
//
{
  double *bg;
  int i;
  int j;
  int k;
  int p;
  double x;
  double y;
  double z;

  bg = new double[3*ng];

  p = 0;

  for ( i = 0; i <= n; i++ )
  {
    x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );
    for ( j = 0; j <= n; j++ )
    {    
      y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );
      for ( k = 0; k <= n; k++ )
      {
        z = c[2] + r * ( double ) ( 2 * k ) / ( double ) ( 2 * n + 1 );

        if ( r * r < pow ( x - c[0], 2 ) 
                   + pow ( y - c[1], 2 )
                   + pow ( z - c[2], 2 ) )
        {
          break;
        }


        bg[0+p*3] = x;
        bg[1+p*3] = y;
        bg[2+p*3] = z;
        p = p + 1;

        if ( 0 < i )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < j )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < k )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < i && 0 < j )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = z;
          p = p + 1;
        }

        if ( 0 < i && 0 < k )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < j && 0 < k )
        {
          bg[0+p*3] = x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }

        if ( 0 < i && 0 < j && 0 < k )
        {
          bg[0+p*3] = 2.0 * c[0] - x;
          bg[1+p*3] = 2.0 * c[1] - y;
          bg[2+p*3] = 2.0 * c[2] - z;
          p = p + 1;
        }
      }
    }
  }

  return bg;
}
//****************************************************************************80

int ball_grid_count ( int n, double r, double c[3] )

//****************************************************************************80
//
//  Purpose:
//
//    BALL_GRID computes grid points inside a ball.
//
//  Discussion:
//
//    The grid is defined by specifying the radius and center of the ball,
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
//    11 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R, the radius of the ball.
//
//    Input, double C[3], the coordinates of the center of the ball.
//
//    Output, int BALL_GRID_COUNT, the number of grid points inside the ball.
//
{
  int i;
  int j;
  int k;
  int ng;
  double x;
  double y;
  double z;

  ng = 0;

  for ( i = 0; i <= n; i++ )
  {
    x = c[0] + r * ( double ) ( 2 * i ) / ( double ) ( 2 * n + 1 );
    for ( j = 0; j <= n; j++ )
    {    
      y = c[1] + r * ( double ) ( 2 * j ) / ( double ) ( 2 * n + 1 );
      for ( k = 0; k <= n; k++ )
      {
        z = c[2] + r * ( double ) ( 2 * k ) / ( double ) ( 2 * n + 1 );

        if ( r * r < pow ( x - c[0], 2 ) 
                   + pow ( y - c[1], 2 )
                   + pow ( z - c[2], 2 ) )
        {
          break;
        }

        ng = ng + 1;

        if ( 0 < i )
        {
          ng = ng + 1;
        }

        if ( 0 < j )
        {
          ng = ng + 1;
        }

        if ( 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < j )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < j && 0 < k )
        {
          ng = ng + 1;
        }

        if ( 0 < i && 0 < j && 0 < k )
        {
          ng = ng + 1;
        }

      }
    }
  }

  return ng;
}
//****************************************************************************80

void r83vec_print_part ( int n, double a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R83VEC_PRINT_PART prints "part" of an R83VEC.
//
//  Discussion:
//
//    An R83VEC is an array of triples of R8's.
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
//    Input, double A[3*N], the vector to be printed.
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
           << "  " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3]  << "\n";
    }
    cout << "  ........  ..............  ..............  ..............\n";
    i = n - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*3]
         << "  " << setw(14) << a[1+i*3] 
         << "  " << setw(14) << a[2+i*3]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      cout << "  " << setw(8) << i
           << ": " << setw(14) << a[0+i*3]
           << "  " << setw(14) << a[1+i*3] 
           << "  " << setw(14) << a[2+i*3]  << "\n";
    }
    i = max_print - 1;
    cout << "  " << setw(8) << i
         << ": " << setw(14) << a[0+i*3]
         << "  " << setw(14) << a[1+i*3] 
         << "  " << setw(14) << a[2+i*3] 
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
