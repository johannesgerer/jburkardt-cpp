# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "line_grid.hpp"

//****************************************************************************80

double *line_grid ( int n, double a, double b, int c )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_GRID: grid points over the interior of a line segment in 1D.
//
//  Discussion:
//
//    In 1D, a grid is created using N points.
//
//    Over the interval [A,B], we have 5 choices for grid centering:
//      1: 0,   1/3, 2/3, 1
//      2: 1/5, 2/5, 3/5, 4/5
//      3: 0,   1/4, 2/4, 3/4
//      4: 1/4, 2/4, 3/4, 1
//      5: 1/8, 3/8, 5/8, 7/8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double A, B, the endpoints.
//
//    Input, int C, the grid centering.
//    1 <= C <= 5.
//
//    Output, double LINE_GRID[N], the points.
//
{
  int j;
  double *x;

  x = new double[n];
//
//  Create the 1D grids in each dimension.
//
  for ( j = 0; j < n; j++ )
  {
    if ( c == 1 )
    {
      if ( n == 1 )
      {
        x[j] = 0.5 * ( a + b );
      }
      else
      {
        x[j] = (   ( double ) ( n - j - 1 ) * a   
                 + ( double ) (     j     ) * b ) 
                 / ( double ) ( n    - 1 );
      }
    }
    else if ( c == 2 )
    {
      x[j] = (   ( double ) ( n - j     ) * a   
               + ( double ) (     j + 1 ) * b ) 
               / ( double ) ( n     + 1 );
    }
    else if ( c == 3 )
    {
      x[j] = (   ( double ) ( n - j     ) * a   
               + ( double ) (     j - 2 ) * b ) 
               / ( double ) ( n         );
    }
    else if ( c == 4 )
    {
      x[j] = (   ( double ) ( n - j - 1 ) * a   
               + ( double ) (     j + 1 ) * b ) 
               / ( double ) ( n         );
    }
    else if ( c == 5 )
    {
      x[j] = (   ( double ) ( 2 * n - 2 * j - 1 ) * a   
               + ( double ) (         2 * j + 1 ) * b ) 
               / ( double ) ( 2 * n             );
    }
  }

  return x;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
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
