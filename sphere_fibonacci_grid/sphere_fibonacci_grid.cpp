# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

# include "sphere_fibonacci_grid.hpp"

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
//    09 November 2014
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
      output << "  " << table[i+j*m];
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//*****************************************************************************/

void sphere_fibonacci_grid_display ( int ng, double xg[], string prefix )

//*****************************************************************************/
//
//  Purpose:
//
//    SPHERE_FIBONACCI_GRID_DISPLAY displays sphere points on a Fibonacci spiral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NG, the number of points.
//
//    Input, double XG[3*NG], the Fibonacci spiral points.
//
//    Input, string PREFIX, a prefix for the filenames.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  string plot_filename;
//
//  Create graphics data file.
//
  data_filename = prefix + "_data.txt";

  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < ng; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      data_unit << "  " << xg[i+j*3];
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created data file '" << data_filename << "'\n";
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
  command_unit << "splot '" << data_filename << "'\n";
  command_unit << "quit\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "%s'\n";

  return;
}
//****************************************************************************80

double *sphere_fibonacci_grid_points ( int ng )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_FIBONACCI_GRID_POINTS computes sphere points on a Fibonacci spiral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 May 2015
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
//    Input, int NG, the number of points.
//
//    Output, double SPHERE_FIBONACCI_GRID_POINTS[3*NG], the Fibonacci 
//    spiral points.
//
{
  double cphi;
  int i;
  double i_r8;
  int j;
  double ng_r8;
  double r8_phi;
  const double r8_pi = 3.141592653589793;
  double sphi;
  double theta;
  double *xyz;

  xyz = new double[3*ng];

  r8_phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;
  ng_r8 = ( double ) ( ng );

  for ( j = 0; j < ng; j++ )
  {
    i_r8 = ( double ) ( - ng + 1 + 2 * j );
    theta = 2.0 * r8_pi * i_r8 / r8_phi;
    sphi = i_r8 / ng_r8;
    cphi = sqrt ( ( ng_r8 + i_r8 ) * ( ng_r8 - i_r8 ) ) / ng_r8;
    xyz[0+j*3] = cphi * sin ( theta );
    xyz[1+j*3] = cphi * cos ( theta );
    xyz[2+j*3] = sphi;
  }

  return xyz;
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
