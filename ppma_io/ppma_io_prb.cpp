# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "ppma_io.hpp"

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
//    PPMA_IO_PRB calls the PPMA_IO test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  timestamp ( );
  cout << "\n";
  cout << "PPMA_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PPMA_IO library.\n";

  error = test01 ( );

  if ( error ) 
  {
    cout << "\n";
    cout << "PPMA_IO_PRB - Fatal error!\n";
    cout << "  TEST01 terminated with an error.\n";
    return 1;
  }

  error = test02 ( );

  if ( error ) 
  {
    cout << "\n";
    cout << "PPMA_IO_PRB - Fatal error!\n";
    cout << "  TEST02 terminated with an error.\n";
    return 1;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PPMA_IO_PRB:\n";
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
//    TEST01 tests PPMA_EXAMPLE and PPMA_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *b;
  bool error;
  string file_out_name = "test01.ascii.ppm";
  int *g;
  int *r;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PPMA_EXAMPLE sets up PPM data.\n";
  cout << "  PPMA_WRITE writes an ASCII PPM file.\n";
  cout << "\n";
  cout << "  Writing the file \"" << file_out_name << "\".\n";

  r = new int[xsize*ysize];
  g = new int[xsize*ysize];
  b = new int[xsize*ysize];

  error = ppma_example ( xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PPMA_EXAMPLE failed!\n";
    return error;
  }

  cout << "\n";
  cout << "  PPMA_EXAMPLE has set up the data.\n";

  error = ppma_write ( file_out_name, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PPMA_WRITE failed!\n";
    return error;
  }

  cout << "\n";
  cout <<  "  PPMA_WRITE was successful.\n";

  delete [] r;
  delete [] g;
  delete [] b;
//
//  Now have PPMA_READ_TEST look at the file we think we created.
//
  error = ppma_read_test ( file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PPMA_READ_TEST failed to read the file we wrote!\n";
    return error;
  }

  cout << "\n";
  cout <<  "  PPMA_READ_TEST was able to read the file we wrote.\n";

  return error;
}
//****************************************************************************80

bool test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PPMA_READ.
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
  string file_in_name = "test02.ascii.ppm";
  int *g;
  int i;
  int j;
  int k;
  int *r;
  int rgb_max;
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PPMA_READ reads the header and data of a PPMA file.\n";
  cout << "\n";
  cout << "  Reading the file \"" << file_in_name << "\".\n";
//
//  Create a data file to read.
//
  error = ppma_write_test ( file_in_name );

  if ( error )
  {
    cout << "\n";
    cout << "TEST02\n";
    cout << "  PPMA_WRITE_TEST failed!\n";
    return error;
  }

  cout << "\n";
  cout << "  PPMA_WRITE_TEST created some test data.\n";
//
//  Now have PPMA_READ try to read it.
//
  error = ppma_read ( file_in_name, xsize, ysize, rgb_max, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "TEST02\n";
    cout << "  PPMA_READ failed!\n";
  }

  cout << "\n";
  cout << "  PPMA_READ read the test data successfully.\n";
  cout << "\n";
  cout << "  Ten sample values:\n";
  cout << "\n";
  for ( k = 0; k < 10; k++ )
  {
    i = ( ( 9 - k ) * 0 + k * ( xsize - 1 ) ) / 9;
    j = ( ( 9 - k ) * 0 + k * ( ysize - 1 ) ) / 9;
    cout << setw(4) << i            << "  "
         << setw(4) << j            << "  "
         << setw(6) << r[i*ysize+j] << "  "
         << setw(6) << g[i*ysize+j] << "  "
         << setw(6) << g[i*ysize+j] << "\n";
  }

  delete [] r;
  delete [] g;
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
