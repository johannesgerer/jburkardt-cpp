# ifdef ANSI_HEADERS
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

int main ( int argc, char *argv[] );
void latinize_test01 ( char *input_filename );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    FSU_LATINIZE_PROBLEMS tests the LATINIZE routines.
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
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "FSU_LATINIZE_PROBLEMS\n";
  cout << "  Test the LATINIZE library.\n";
  cout << "\n";
  cout << "  Read a dataset of N points in M dimensions,\n";
  cout << "  modify it into a Latin hypercube,\n";
  cout << "  write the modified dataset to a file.\n";

  latinize_test01 ( "cvt_07_00010.txt" );
  latinize_test01 ( "halton_02_00100.txt" );

  cout << "\n";
  cout << "FSU_LATINIZE_PROBLEMS\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void latinize_test01 ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    LATINIZE_TEST01 latinizes the data in a given file.
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
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  ofstream output;
  char *output_filename;
  double *table;
//
//  Need to create the output file name from the input filename.
//
  output_filename = file_name_ext_swap ( input_filename, "latin.txt" );

  dtable_header_read ( input_filename, &m, &n );

  cout << "\n";
  cout << "  Read the header of \"" << input_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Number of points N  = " << n << "\n";

  table = dtable_data_read ( input_filename, m, n );

  cout << "\n";
  cout << "  Read the data in \"" << input_filename << "\".\n";

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of data read from file:" );

  fsu_latinize ( m, n, table );

  cout << "\n";
  cout << "  Latinized the data.\n";

  r8mat_transpose_print_some ( m, n, table, 1, 1, 5, 5, 
    "  Small portion of Latinized data:" );
//
//  Write the data to a file.
//
  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "LATINIZE_TEST01 - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  dtable_header_write ( m, n, output_filename, output );

  output << "#  The data was read from " << input_filename 
         << " by DTABLE_DATA_READ.CC\n";
  output << "#  The data was latinized by FSU_LATINIZE.C.\n";
  output << "#  The data was written to " << output_filename 
         << " by DTABLE_DATA_WRITE.CC.\n";
  output << "#\n";

  dtable_data_write ( m, n, table, output );

  output.close ( );

  cout << "\n";
  cout << "  Wrote the latinized data to \"" << output_filename << "\".\n";

  delete [] table;

  return;
}
