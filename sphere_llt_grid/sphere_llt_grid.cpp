# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

# include "sphere_llt_grid.hpp"

//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
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
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
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
    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j - 1 << ":";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

void sphere_llt_grid_display ( int ng, double xg[], int line_num, 
  int line_data[], string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_DISPLAY displays a latitude/longitude triangle grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NG, the number of points.
//
//    Input, double XG[3*NG], the points.
//
//    Input, int LINE_NUM, the number of grid lines.
//
//    Input, inte LINE_DATA[2*LINE_NUM], contains pairs of 
//    point indices for line segments that make up the grid.
//
//    Input, string PREFIX, a prefix for the filenames.
//
{
  string command_filename;
  ofstream command_unit;
  int i;
  int j;
  int j1;
  int j2;
  int l;
  string line_filename;
  ofstream line_unit;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
//
//  Create graphics data files.
//
  node_filename = prefix + "_nodes.txt";

  node_unit.open ( node_filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      node_unit << "  " << xg[i+j*3];
    }
    node_unit <<"\n";
  }
  node_unit.close ( );
  cout << "\n";
  cout << "  Created node file '" << node_filename << "'\n";

  line_filename = prefix + "_lines.txt";

  line_unit.open ( line_filename.c_str ( ) );
  for ( l = 0; l < line_num; l++ )
  {
    if ( 0 < l )
    {
      line_unit << "\n";
      line_unit << "\n";
    }
    j1 = line_data[0+l*2];
    j2 = line_data[1+l*2];
    line_unit << "  " << xg[0+j1*3]
              << "  " << xg[1+j1*3]
              << "  " << xg[2+j1*3] << "\n";
    line_unit << "  " << xg[0+j2*3]
              << "  " << xg[1+j2*3]
              << "  " << xg[2+j2*3] << "\n";
  }
  line_unit.close ( );
  cout << "\n";
  cout << "  Created line file '" << line_filename << "'\n";
//
//  Create graphics command file.
//
  command_filename = prefix + "_commands.txt";

  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";

  plot_filename = prefix + ".png";

  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << prefix << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set style data points\n";
  command_unit << "set timestamp\n";
  command_unit << "set view equal xyz\n";
  command_unit << "splot '" << line_filename << "' with lines lw 3, \\\n";
  command_unit << "      '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit << "quit\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

int sphere_llt_grid_line_count ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_LINE_COUNT counts latitude/longitude triangle grid lines.
//
//  Discussion:
//
//    A SPHERE LLT grid imposes a grid of triangles on a sphere,
//    using latitude and longitude lines.
//
//    The number returned is the number of pairs of points to be connected.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude and
//    longitude lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so LAT_NUM = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_LLT_GRID_LINE_COUNT, the number of grid lines.
//
{
  int line_num;

  line_num = long_num * ( lat_num + 1 ) 
           + long_num * lat_num
           + long_num * ( lat_num - 1 );

  return line_num;
}
//****************************************************************************80

int *sphere_llt_grid_lines ( int nlat, int nlong, int line_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_LINES: latitude/longitude triangle grid lines.
//
//  Discussion:
//
//    A SPHERE LLT grid imposes a grid of triangles on a sphere,
//    using latitude and longitude lines.
//
//    The point numbering system is the same used in SPHERE_LLT_POINTS,
//    and that routine may be used to compute the coordinates of the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NLAT, NLONG, the number of latitude and longitude
//    lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so NLAT = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Input, int LINE_NUM, the number of grid lines.
//
//    Output, int SPHERE_LLT_GRID_LINES[2*LINE_NUM], contains pairs of point 
//    indices for line segments that make up the grid.
//
{
  int i;
  int j;
  int l;
  int *line;
  int next;
  int newcol;
  int old;

  line = new int[2*line_num];
  l = 0;
//
//  "Vertical" lines.
//
  for ( j = 0; j <= nlong - 1; j++ )
  {
    old = 0;
    next = j + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;

    for ( i = 1; i <= nlat-1; i++ )
    {
      old = next;
      next = old + nlong;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    line[0+l*2] = old;
    line[1+l*2] = 1 + nlat * nlong;
    l = l + 1;
  }
//
//  "Horizontal" lines.
//
  for ( i = 1; i <= nlat; i++ )
  {
    next = ( i - 1 ) * nlong + 1;

    for ( j = 0; j <= nlong - 2; j++ )
    {
      old = next;
      next = old + 1;
      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }

    old = next;
    next = ( i - 1 ) * nlong + 1;
    line[0+l*2] = old;
    line[1+l*2] = next;
    l = l + 1;
  }
//
//  "Diagonal" lines.
//
  for ( j = 0; j < nlong; j++ )
  {
    old = 0;
    next = j + 1;
    newcol = j;

    for ( i = 1; i < nlat; i++ )
    {
      old = next;
      next = old + nlong + 1;
      newcol = newcol + 1;
      if ( nlong - 1 < newcol )
      {
        newcol = 0;
        next = next - nlong;
      }

      line[0+l*2] = old;
      line[1+l*2] = next;
      l = l + 1;
    }
  }
  return line;
}
//****************************************************************************80

int sphere_llt_grid_point_count ( int lat_num, int long_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_POINT_COUNT counts points for a latitude/longitude grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int LAT_NUM, LONG_NUM, the number of latitude 
//    and longitude lines to draw.  The latitudes do not include the North and 
//    South poles, which will be included automatically, so LAT_NUM = 5, for 
//    instance, will result in points along 7 lines of latitude.
//
//    Output, int SPHERE_LLT_GRID_POINT_COUNT, the number of grid points.
//
{
  int point_num;

  point_num = 2 + lat_num * long_num;

  return point_num;
}
//****************************************************************************80

double *sphere_llt_grid_points ( double r, double pc[3], int lat_num, 
  int lon_num, int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_LLT_GRID_POINTS produces points on a latitude/longitude grid.
//
//  Discussion:
//
//    A SPHERE LLT grid imposes a grid of triangles on a sphere,
//    using latitude and longitude lines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the radius of the sphere.
//
//    Input, double PC[3], the coordinates of the center of the sphere.
//
//    Input, int LAT_NUM, LON_NUM, the number of latitude and longitude
//    lines to draw.  The latitudes do not include the North and South
//    poles, which will be included automatically, so LAT_NUM = 5, for instance,
//    will result in points along 7 lines of latitude.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double SPHERE_LLT_GRID_POINTS[3*POINT_NUM], the coordinates 
//    of the grid points.
//
{
  int lat;
  int lon;
  int n;
  double *p;
  double phi;
  double r8_pi = 3.141592653589793;
  double theta;

  p = new double[3*point_num];
  n = 0;
//
//  The north pole.
//
  theta = 0.0;
  phi = 0.0;

  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;
//
//  Do each intermediate ring of latitude.
//
  for ( lat = 1; lat <= lat_num; lat++ )
  {
    phi = ( double ) ( lat ) * r8_pi / ( double ) ( lat_num + 1 );
//
//  Along that ring of latitude, compute points at various longitudes.
//
    for ( lon = 0; lon < lon_num; lon++ )
    {
      theta = ( double ) ( lon ) * 2.0 * r8_pi / ( double ) ( lon_num );

      p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
      p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
      p[2+n*3] = pc[2] + r * cos ( phi );
      n = n + 1;
    }
  }
//
//  The south pole.
//
  theta = 0.0;
  phi = r8_pi;
  p[0+n*3] = pc[0] + r * sin ( phi ) * cos ( theta );
  p[1+n*3] = pc[1] + r * sin ( phi ) * sin ( theta );
  p[2+n*3] = pc[2] + r * cos ( phi );
  n = n + 1;

  return p;
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
