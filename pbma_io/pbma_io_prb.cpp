# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "pbma_io.hpp"

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PBMA_IO_PRB.
//
//  Discussion:
//
//    PBMA_IO_PRB calls the PBMA_IO test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PBMA_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PBMA_IO library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PBMA_IO_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PBMA_EXAMPLE, PBMA_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  string file_out_name = "pbma_io_prb_01.ascii.pbm";
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PBMA_EXAMPLE sets up ASCII PBM data.\n";
  cout << "  PBMA_WRITE writes an ASCII PBM file.\n";
  cout << "\n";
  cout << "  Writing the file \"" << file_out_name << "\".\n";

  b = new int[xsize*ysize];

  pbma_example ( xsize, ysize, b );

  cout << "\n";
  cout << "  PBMA_EXAMPLE has set up the data.\n";

  pbma_write ( file_out_name, xsize, ysize, b );

  cout << "\n";
  cout <<  "  PBMA_WRITE was successful.\n";

  delete [] b;
//
//  Now have PBMA_READ_TEST look at the file we think we created.
//
  pbma_read_test ( file_out_name );

  cout << "\n";
  cout <<  "  PBMA_READ_TEST was able to read the file.\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PBMA_READ_HEADER, PBMA_READ_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  string file_in_name = "pbma_io_prb_02.ascii.pbm";
  int i;
  int j;
  int k;
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PBMA_READ reads the header and data of an ASCII PBM file.\n";
//
//  Create a data file to read.
//
  pbma_write_test ( file_in_name );

  cout << "\n";
  cout << "  PBMA_WRITE_TEST created some test data.\n";
//
//  Now have PBMA_READ try to read it.
//
  pbma_read ( file_in_name, xsize, ysize, &b );

  cout << "\n";
  cout << "  PBMA_READ was able to read the file we created.\n";
  cout << "\n";
  cout << "  Sample data:\n";
  cout << "\n";

  for ( k = 0; k <= 29; k++ )
  {
    i = ( ( 29 - k ) * 0 + k * ( xsize - 1 ) ) / 29;
    j = ( ( 29 - k ) * 0 + k * ( ysize - 1 ) ) / 29;
    cout << setw(4) << i            << "  "
         << setw(4) << j            << "  "
         << setw(6) << b[i*ysize+j] << "\n";
  }

  delete [] b;

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
