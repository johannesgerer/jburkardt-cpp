# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "pgmb_io.hpp"

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
//    MAIN tests the PGMB_IO_PRB library.
//
//  Discussion:
//
//    PGMB_IO_PRB tests the PGMB_IO library.
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
  cout << "PGMB_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PGMB_IO library.\n";

  error = test01 ( );

  if ( error ) 
  {
    cout << "\n";
    cout << "PGMB_IO_PRB - Fatal error!\n";
    cout << "  TEST01 terminated with an error.\n";
    return 1;
  }

  error = test02 ( );

  if ( error ) 
  {
    cout << "\n";
    cout << "PGMB_IO_PRB - Fatal error!\n";
    cout << "  TEST02 terminated with an error.\n";
    return 1;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PGMB_IO_PRB:\n";
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
//    TEST01 tests PGMB_EXAMPLE, PGMB_WRITE.
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
  string file_out_name = "pgmb_io_prb_01.pgm";
  unsigned char *g;
  int i;
  unsigned char *indexg;
  int j;
  unsigned char maxg;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PGMB_EXAMPLE sets up PGMB data.\n";
  cout << "  PGMB_WRITE writes a PGMB file.\n";
  cout << "\n";
  cout << "  Writing the file \"" << file_out_name << "\".\n";

  g = new unsigned char[xsize*ysize];

  error = pgmb_example ( xsize, ysize, g );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PGMB_EXAMPLE failed!\n";
    return error;
  }
  else
  {
    cout << "\n";
    cout << "  PGMB_EXAMPLE has set up the data.\n";
  }

  maxg = 0;
  indexg = g;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      if ( maxg < *indexg )
      {
        maxg = *indexg;
      }
      indexg = indexg + 1;
    }
  }
  cout << "\n";
  cout << "  Gray scale data has maximum value " << ( int ) maxg << "\n";

  error = pgmb_write ( file_out_name, xsize, ysize, g );

  if ( error )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  PGMB_WRITE failed!\n";
  }
  else
  {
    cout << "\n";
    cout <<  "  PGMB_WRITE was successful.\n";
  }

  delete [] g;
//
//  Now have PGMB_READ_TEST look at the file we think we created.
//
  error = pgmb_read_test ( file_out_name );

  return error;
}
//****************************************************************************80

bool test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PGMB_READ_HEADER, PGMB_READ_DATA.
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
  string file_in_name = "pgmb_io_prb_02.pgm";
  unsigned char *g;
  unsigned char maxg;
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PGMB_READ reads the header and data of a PGMB file.\n";
  cout << "\n";
  cout << "  Reading the file \"" << file_in_name << "\".\n";
//
//  Create a data file to read.
//
  error = pgmb_write_test ( file_in_name );

  if ( error )
  {
    cout << "\n";
    cout << "  PGMB_WRITE_TEST failed!\n";
  }
  else
  {
    cout << "\n";
    cout << "  PGMB_WRITE_TEST created some test data.\n";
  }
//
//  Now have PGMB_READ try to read it.
//
  error = pgmb_read ( file_in_name, xsize, ysize, maxg, &g );

  if ( error )
  {
    cout << "\n";
    cout << "  PGMB_READ failed!\n";
  }
  else
  {
    cout << "\n";
    cout << "  PGMB_READ read the test data successfully.\n";
  }

  delete [] g;

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
