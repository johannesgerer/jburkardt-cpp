# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "pbmb_io.hpp"

int main ( int argc, char *argv[] );
bool test01 ( );
bool test02 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PBMB_IO_PRB.
//
//  Discussion:
//
//    PBMB_IO_PRB calls the PBMB_IO test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  timestamp ( );
  cout << "\n";
  cout << "PBMB_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PBMB_IO library.\n";

  error = test01 ( );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_IO_PRB - Fatal error!\n";
    cout << "  TEST01 terminated with an error.\n";
    return 1;
  }

  error = test02 ( );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_IO_PRB - Fatal error!\n";
    cout << "  TEST02 terminated with an error.\n";
    return 1;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PBMB_IO_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PBMB_EXAMPLE, PBMB_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  bool error;
  string file_out_name = "pbmb_io_prb_01.pbm";
  unsigned char *indexb;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PBMB_EXAMPLE sets up PBMB data.\n";
  cout << "  PBMB_WRITE writes a PBMB file.\n";
  cout << "\n";
  cout << "  Writing the file \"" << file_out_name << "\".\n";

  b = new int[xsize*ysize];

  error = pbmb_example ( xsize, ysize, b );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PBMB_EXAMPLE failed!\n";
    return error;
  }
  else
  {
    cout << "\n";
    cout << "  PBMB_EXAMPLE has set up the data.\n";
  }

  error = pbmb_write ( file_out_name, xsize, ysize, b );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PBMB_WRITE failed!\n";
  }
  else
  {
    cout << "\n";
    cout <<  "  PBMB_WRITE was successful.\n";
  }

  delete [] b;
//
//  Now have PBMB_READ_TEST look at the file we think we created.
//
  error = pbmb_read_test ( file_out_name );

  return error;
}
//****************************************************************************80

bool test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PBMB_READ_HEADER, PBMB_READ_DATA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  bool error;
  string file_in_name = "pbmb_io_prb_02.pbm";
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PBMB_READ reads the header and data of a PBMB file.\n";
  cout << "\n";
  cout << "  Reading the file \"" << file_in_name << "\".\n";
//
//  Create a data file to read.
//
  error = pbmb_write_test ( file_in_name );

  if ( error )
  {
    cout << "\n";
    cout << "  PBMB_WRITE_TEST failed!\n";
  }
  else
  {
    cout << "\n";
    cout << "  PBMB_WRITE_TEST created some test data.\n";
  }
//
//  Now have PBMB_READ try to read it.
//
  error = pbmb_read ( file_in_name, xsize, ysize, &b );

  if ( error )
  {
    cout << "\n";
    cout << "  PBMB_READ failed!\n";
  }
  else
  {
    cout << "\n";
    cout << "  PBMB_READ read the test data successfully.\n";
  }

  delete [] b;

  return error;
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
