# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "sparse_display.hpp"

//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
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

void spy_file ( string header, string data_filename )

//****************************************************************************80
//
//  Purpose:
//
//    SPY_FILE plots a sparsity pattern stored in a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string HEADER, the name to be used for the
//    title of the plot, and as part of the names of the command
//    and plot files.
//
//    Input, string DATA_FILENAME, the name of the file
//    containing the indices of nonzero matrix entries.
//
{
  string command_filename;
  ofstream command_unit;
  ifstream data_unit;
  int i;
  const int i4_huge = 2147483647;
  int j;
  int m0;
  int m1;
  int n0;
  int n1;
  int nz_num;
  string png_filename;

  n0 = + i4_huge;
  n1 = - i4_huge;
  m0 = + i4_huge;
  m1 = - i4_huge;
  nz_num = 0;

  data_unit.open ( data_filename.c_str ( ) );

  for ( ; ; )
  {
    data_unit >> i;
    if ( data_unit.eof ( ) )
    {
      break;
    }

    data_unit >> j;
    if ( data_unit.eof ( ) )
    {
      break;
    }

    nz_num = nz_num + 1;
    m0 = i4_min ( m0, i );
    m1 = i4_max ( m1, i );
    n0 = i4_min ( n0, j );
    n1 = i4_max ( n1, j );
  }

  data_unit.close ( );
//
//  Create command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "unset key\n";
  command_unit << "set term png\n";

  png_filename = header + ".png";
  command_unit << "set output '" << png_filename << "'\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set xlabel '<--- J --->'\n";
  command_unit << "set ylabel '<--- I --->'\n";
  
  command_unit << "set title '" 
               << nz_num << " nonzeros for \""
               << header << "\"'\n";
  command_unit << "set timestamp\n";
  command_unit << "plot [y="
               << m0 << ":"
               << m1 << "] [x="
               << n0 << ":"
               << n1 << "] '"
               << data_filename << "' with points pt 5\n";
 
  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

void spy_ge ( int m, int n, double a[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    SPY_GE plots a sparsity pattern for a general storage (GE) matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns
//    in the matrix.
//
//    Input, double A[M*N], the matrix.
//
//    Input, string HEADER, the name to be used for the
//    title of the plot, and as part of the names of the data, command
//    and plot files.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  int nz_num;
  string png_filename;
//
//  Create data file.
//
  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  nz_num = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] != 0.0 )
      {
        data_unit << j << "  "
                  << i << "\n";
        nz_num = nz_num + 1;
      }
    }
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created sparsity data file '" << data_filename << "'\n";
//
//  Create command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "unset key\n";
  command_unit << "set term png\n";

  png_filename = header + ".png";
  command_unit << "set output '" << png_filename << "'\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set xlabel '<--- J --->'\n";
  command_unit << "set ylabel '<--- I --->'\n";
  command_unit << "set title '" << nz_num << " nonzeros for \"" 
               << header << "\"'\n";
  command_unit << "set timestamp\n";
  command_unit << "plot [y=0:" << n - 1 << "] [x="
               << m - 1 << ":0] '"
               << data_filename << "' with points pt 5\n";
 
  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

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
