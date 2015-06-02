# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

# include "latinize.hpp"

int main ( int argc, char *argv[] );
void test01 ( string input_filename );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LATINIZE_PRB.
//
//  Discussion:
//
//    LATINIZE_PRB tests the LATINIZE library.
//
//    The dataset is presumed to be an M by N array of real numbers,
//    where M is the spatial dimension, and N is the number of sample points.
//
//    The dataset is presumed to be stored in a file, with N records,
//    one per each sample point.  (Comment records may be included, 
//    which begin with '#'.)
//
//    The program reads the data file, "latinizes" the data, and writes
//    the latinized data to a new file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2004
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LATINIZE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LATINIZE library.\n";
  cout << "\n";
  cout << "  Read a dataset of N points in M dimensions,\n";
  cout << "  modify it into a Latin hypercube,\n";
  cout << "  write the modified dataset to a file.\n";

  test01 ( "cvt_02_00010.txt" );
  test01 ( "cvt_03_00007.txt" );
  test01 ( "cvt_03_00056.txt" );
  test01 ( "cvt_07_00100.txt" );
//
//  Terminate.
//
  cout << "\n";
  cout << "LATINIZE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    LATINIZE_PRB tests the LATINIZE routines.
//
//  Discussion:
//
//    The dataset is presumed to be an M by N array of real numbers,
//    where M is the spatial dimension, and N is the number of sample points.
//
//    The dataset is presumed to be stored in a file, with N records,
//    one per each sample point.  (Comment records may be included, 
//    which begin with '#'.)
//
//    The program reads the data file, "latinizes" the data, and writes
//    the latinized data to a new file.
//
//  Modified:
//
//    08 October 2004
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  ofstream output;
  string output_filename;
  double *table;
//
//  Need to create the output file name from the input filename.
//
  output_filename = file_name_ext_swap ( input_filename, "latin.txt" );

  r8mat_header_read ( input_filename, &m, &n );

  cout << "\n";
  cout << "  Read the header of \"" << input_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of points N  = " << n << "\n";

  table = r8mat_data_read ( input_filename, m, n );

  cout << "\n";
  cout << "  Read the data in \"" << input_filename << "\".\n";

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of data read from file:" );

  r8mat_latinize ( m, n, table );

  cout << "\n";
  cout << "  Latinized the data.\n";

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of Latinized data:" );
//
//  Write the data to a file.
//
  r8mat_write ( output_filename, m, n, table );

  cout << "\n";
  cout << "  Wrote the latinized data to \"" << output_filename << "\".\n";

  delete [] table;

  return;
}
