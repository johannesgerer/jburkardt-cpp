# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iostream>
# include <iomanip>

using namespace std;

# include "polygon_grid.hpp"

//****************************************************************************80

int polygon_grid_count ( int n, int nv )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_COUNT counts the grid points inside a polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals on a side.
//
//    Input, int NV, the number of vertices.
//    3 <= NV.
//
//    Output, int POLYGON_GRID_COUNT, the number of grid points.
//
{
  int ng;

  ng = 1 + nv * ( n * ( n + 1 ) ) / 2;

  return ng;
}
//****************************************************************************80

void polygon_grid_display ( int n, int nv, double v[], int ng, double xg[], 
  string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_DISPLAY displays grid points inside a polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NV, the number of vertices in the polygon.
//
//    Input, double V[2*NV], the coordinates of the vertices.
//
//    Input, int NG, the number of grid points.
//
//    Input, double XG[2*NG], the grid points.
//
//    Input, string PREFIX, a string used to name the files.
//
{
  ofstream command_unit;
  string command_filename;
  ofstream grid_unit;
  string grid_filename;
  int j;
  string plot_filename;
  double vc[2];
  ofstream vertex_unit;
  string vertex_filename;
//
//  Determine the centroid.
//
  vc[0] = 0.0;
  vc[1] = 0.0;
  for ( j = 0; j < nv; j++ )
  {
    vc[0] = vc[0] + v[0+j*2];
    vc[1] = vc[1] + v[1+j*2];
  }
  vc[0] = vc[0] / ( double ) ( nv );
  vc[1] = vc[1] / ( double ) ( nv );
//
//  Write the vertex file.
//
  vertex_filename = prefix + "_vertex.txt";
  vertex_unit.open ( vertex_filename.c_str ( ) );

  for ( j = 0; j < nv; j++ )
  {
    vertex_unit << "  " << v[0+2*j]
                << "  " << v[1+2*j] << "\n";
  }
  vertex_unit << "  " << v[0+0*j]
              << "  " << v[1+0*j] << "\n";
  for ( j = 0; j < nv; j++ )
  {
    vertex_unit << "\n";
    vertex_unit << "  " << v[0+j*2]
                << "  " << v[1+j*2] << "\n";
    vertex_unit << "  " << vc[0]
                << "  " << vc[1] << "\n";
  }
  vertex_unit.close ( );
  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Write the gridpoint file.
//
  grid_filename = prefix + "_grid.txt";
  grid_unit.open ( grid_filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    grid_unit << "  " << xg[0+j*2]
              << "  " << xg[1+j*2] << "\n";
  }
  grid_unit.close ( );
  cout << "\n";
  cout << "  Created grid file '" << grid_filename << "'\n";
//
//  Write the command file.
//
  plot_filename = prefix + ".png";

  command_filename = prefix + "_commands.txt";

  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title '" << prefix << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set style data lines\n";

  command_unit << 
    "plot '" << grid_filename << "' using 1:2 with points lt 3 pt 3,\\\n";
  command_unit << 
    "     '" << vertex_filename << "' using 1:2 lw 3 linecolor rgb 'black'\n";
  command_unit << "quit\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

double *polygon_grid_points ( int n, int nv, double v[], int ng )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_GRID_POINTS computes points on a polygonal grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, int NV, the number of vertices in the polygon.
//
//    Input, double V[2*NV], the coordinates of the vertices.
//
//    Input, int NG, the number of grid points.
//
//    Output, double POLYGON_GRID_POINTS[2*NG], the coordinates of the 
//    grid points.
//
{
  int i;
  int j;
  int k;
  int l;
  int lp1;
  int p;
  double vc[2];
  double *xg;

  xg = new double[2*ng];
  p = 0;
//
//  Determine the centroid.
//
  vc[0] = 0.0;
  vc[1] = 0.0;
  for ( j = 0; j < nv; j++ )
  {
    vc[0] = vc[0] + v[0+j*2];
    vc[1] = vc[1] + v[1+j*2];
  }
  vc[0] = vc[0] / ( double ) ( nv );
  vc[1] = vc[1] / ( double ) ( nv );
//
//  The centroid is the first point.
//
  xg[0+p*2] = vc[0];
  xg[1+p*2] = vc[1];
  p = p + 1;
//
//  Consider each triangle formed by two consecutive vertices and the centroid,
//  but skip the first line of points.
//
  for ( l = 0; l < nv; l++ )
  {
    lp1 = ( ( l + 1 ) % nv );
    for ( i = 1; i <= n; i++ )
    {
      for ( j = 0; j <= n - i; j++ )
      {
        k = n - i - j;
        xg[0+p*2] = ( ( double ) ( i ) * v[0+l*2]   
                    + ( double ) ( j ) * v[0+lp1*2] 
                    + ( double ) ( k ) * vc[0] )  
                    / ( double ) ( n );
        xg[1+p*2] = ( ( double ) ( i ) * v[1+l*2]   
                    + ( double ) ( j ) * v[1+lp1*2] 
                    + ( double ) ( k ) * vc[1] )  
                    / ( double ) ( n );
        p = p + 1;
      }
    }
  }

  return xg;
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
