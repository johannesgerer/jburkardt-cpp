# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>
# include <fstream>

using namespace std;

# include "pyramid_grid.hpp"

//****************************************************************************80

int pyramid_grid_size ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_GRID_SIZE sizes a pyramid grid.
//
//  Discussion:
//
//    0:  x
//
//    1:  x  x
//        x  x
//
//    2:  x  x  x
//        x  x  x
//        x  x  x
//
//    3:  x  x  x  x
//        x  x  x  x
//        x  x  x  x
//        x  x  x  x
//
//    N  Size
//
//    0     1
//    1     5 = 1 + 4
//    2    14 = 1 + 4 + 9
//    3    30 = 1 + 4 + 9 + 16
//    4    55 = 1 + 4 + 9 + 16 + 25
//    5    91 = 1 + 4 + 9 + 16 + 25 + 36
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Output, int PYRAMID_GRID_SIZE, the number of
//    nodes in the grid of size N.
//
{
  int np1;
  int value;

  np1 = n + 1;

  value = ( np1 * ( np1 + 1 ) * ( 2 * np1 + 1 ) ) / 6;

  return value;
}
//****************************************************************************80

double *pyramid_unit_grid ( int n, int ng )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_GRID computes grid points in the unit pyramid.
//
//  Discussion:
//
//    The unit pyramid has base (-1,-1,0), (+1,1,0), (+1,+1,0), (-1,+1,0)
//    and vertex (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NG, the number of nodes to generate,
//    as determined by pyramid_grid_size().
//
//    Output, double PYRAMID_UNIT_GRID[3*NG], the grid point coordinates.
//
{
  int g;
  int hi;
  int i;
  int j;
  int k;
  int lo;
  double *pg;

  pg = new double[3*ng];

  g = 0;

  for ( k = n; 0 <= k; k-- )
  {
    hi = n - k;
    lo = - hi;
    for ( j = lo; j <= hi; j = j + 2 )
    {
      for ( i = lo; i <= hi; i = i + 2 )
      {
        pg[0+g*3] = ( double ) ( i ) / ( double ) ( n );
        pg[1+g*3] = ( double ) ( j ) / ( double ) ( n );
        pg[2+g*3] = ( double ) ( k ) / ( double ) ( n );
        g = g + 1;
      }
    }
  }

  return pg;
}
//****************************************************************************80

void pyramid_unit_grid_plot ( int n, int ng, double pg[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_GRID_PLOT sets up a GNUPLOT plot of a unit pyramid grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NG, the number of nodes to generate,
//    as determined by pyramid_grid_size().
//
//    Input, double PG[3*NG], the grid point coordinates.
//
//    Input, string HEADER, the header for the files.
//
{
  string command_filename;
  ofstream command_unit;
  int j;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
  double v1[3];
  double v2[3];
  double v3[3];
  double v4[3];
  double v5[3];
  string vertex_filename;
  ofstream vertex_unit;
//
//  Create the vertex file.
//
  pyramid_unit_vertices ( v1, v2, v3, v4, v5 );

  vertex_filename = header + "_vertices.txt";

  vertex_unit.open ( vertex_filename.c_str ( ) );

  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << v5[0] << "  "
              << v5[1] << "  "
              << v5[2] << "\n";
  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v2[0] << "  "
              << v2[1] << "  "
              << v2[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v3[0] << "  "
              << v3[1] << "  "
              << v3[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v4[0] << "  "
              << v4[1] << "  "
              << v4[2] << "\n";
  vertex_unit << "\n";

  vertex_unit << v1[0] << "  "
              << v1[1] << "  "
              << v1[2] << "\n";
  vertex_unit << v5[0] << "  "
              << v5[1] << "  "
              << v5[2] << "\n";
  vertex_unit << "\n";

  vertex_unit.close ( );

  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Create the node file.
//
  node_filename = header + "_nodes.txt";

  node_unit.open ( node_filename.c_str ( ) );

  for ( j = 0; j < ng; j++ )
  {
    node_unit << pg[0+j*3] << "  "
              << pg[1+j*3] << "  "
              << pg[2+j*3] << "\n";
  }
  node_unit.close ( );

  cout << " Created node file '" << node_filename << "'\n";
//
//  Create the command file.
//
  command_filename = header + "_commands.txt";

  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";

  plot_filename = header + ".png";

  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set view equal xyz\n";
  command_unit << "set view 80, 40\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "splot '" << vertex_filename << "' with lines lw 3, \\\n";
  command_unit << "      '" << node_filename << "' with points pt 7 lt 0\n";

  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

void pyramid_unit_vertices ( double v1[], double v2[], double v3[], 
  double v4[], double v5[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID_UNIT_VERTICES returns the vertices of the unit pyramid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], the vertices.
//
{
  v1[0] =  0.0;
  v1[1] =  0.0;
  v1[2] = +1.0;

  v2[0] = -1.0;
  v2[1] = -1.0;
  v2[2] =  0.0;

  v3[0] = +1.0;
  v3[1] = -1.0;
  v3[2] =  0.0;

  v4[0] = +1.0;
  v4[1] = +1.0;
  v4[2] =  0.0;

  v5[0] = -1.0;
  v5[1] = +1.0;
  v5[2] =  0.0;

  return;
}
//****************************************************************************80

void r8_print ( double r, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PRINT prints an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the value to print.
//
//    Input, string TITLE, a title.
//
{
  cout << title << "  "
       << r << "\n";

  return;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
