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
void handle ( char *input_filename );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TABLE_LATINIZE.
//
//  Discussion:
//
//    TABLE_LATINIZE is the main routine of a program to "latinize" a dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    table_latinize input_filename
//
{
  int i;
  char input_filename[80];

  timestamp ( );

  cout << "\n";
  cout << "LATINIZE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Read a dataset of N points in M dimensions,\n";
  cout << "  modify it into a Latin hypercube,\n";
  cout << "  write the modified dataset to a file.\n";

  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TABLE_LATINIZE:\n";
    cout << "  Please enter the name of a file to be analyzed.\n";

    cin.getline ( input_filename, sizeof ( input_filename ) );

    handle ( input_filename );

  }
  else 
  {
    for ( i = 1; i < argc; i++ ) 
    {
      handle ( argv[i] );
    }
  } 

  cout << "\n";
  cout << "TABLE_LATINIZE\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void handle ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE handles a single file.
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
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//  Local parameters:
//
//    Local, int DIM_NUM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[DIM_NUM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
{
  int dim_num;
  int n;
  ofstream output;
  char *output_filename;
  double *table;
//
//  Need to create the output file name from the input filename.
//
  output_filename = file_name_ext_swap ( input_filename, "latin.txt" );

  dtable_header_read ( input_filename, &dim_num, &n );

  cout << "\n";
  cout << "  Read the header of \"" << input_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  Number of points N  = " << n << "\n";

  table = dtable_data_read ( input_filename, dim_num, n );

  cout << "\n";
  cout << "  Read the data in \"" << input_filename << "\".\n";

  dtable_print_some ( dim_num, n, table, 1, 1, 5, 5, 
    "  Small portion of data read from file:" );

  dtable_latinize ( dim_num, n, table );

  cout << "\n";
  cout << "  Latinized the data.\n";

  dtable_print_some ( dim_num, n, table, 1, 1, 5, 5, 
    "  Small portion of Latinized data:" );
//
//  Write the data to a file.
//
  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  dtable_header_write ( dim_num, n, output_filename, output );

  output << "#  This data was created by the C++ TABLE_LATINIZE program.\n";
  output << "#  The data was read from " << input_filename 
         << " by DTABLE_DATA_READ.C\n";
  output << "#  The data was latinized by DTABLE_LATINIZE.CC.\n";
  output << "#  The data was written to " << output_filename 
         << " by DTABLE_DATA_WRITE.C.\n";
  output << "#\n";

  dtable_data_write ( dim_num, n, table, output );

  output.close ( );

  cout << "\n";
  cout << "  Wrote the latinized data to \"" << output_filename << "\".\n";

  delete [] table;

  return;
}
