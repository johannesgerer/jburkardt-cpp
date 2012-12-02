# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "ppmb_io.hpp"

int main ( int argc, char *argv[] );
bool test01 ( string file_name );
bool test02 ( string file_name );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN calls the PPMB_IO test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  timestamp ( );
  cout << "\n";
  cout << "PPMB_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PPMB_IO library.\n";

  error = test01 ( "ppmb_io_prb.ppm" );
  error = test02 ( "ppmb_io_prb.ppm" );
//
//  Terminate.
//
  cout << "\n";
  cout << "PPMB_IO_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool test01 ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PPM_EXAMPLE, PPMB_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 September 2008
//
//  Author:
//
//    John Burkardt
//
{
  unsigned char *b;
  bool error;
  unsigned char *g;
  unsigned char *r;
  int result;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PPM_EXAMPLE sets up PPM data.\n";
  cout << "  PPMB_WRITE writes a binary PPM file.\n";

  r = new unsigned char[ xsize * ysize ];
  g = new unsigned char[ xsize * ysize ];
  b = new unsigned char[ xsize * ysize ];

  error = ppmb_example ( xsize, ysize, r, g, b );

  cout << "\n";
  cout << "  PPM_EXAMPLE has set up the data.\n";

  error = ppmb_write ( file_name, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "  PPMB_WRITE failed!\n";
  }
  else
  {
    cout << "  PPMB_WRITE was successful.\n";
  }

  delete [] r;
  delete [] g;
  delete [] b;

  return error;
}
//****************************************************************************80

bool test02 ( string file_in_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PPMB_READ_HEADER, PPMB_READ_DATA.
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
  unsigned char *b;
  bool error;
  ifstream file_in;
  unsigned char *g;
  int maxrgb;
  unsigned char *r;
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PPMB_READ_HEADER reads the header of a PPMB file.\n";
  cout << "  PPMB_READ_DATA reads the data in a PPMB file.\n";
 
  cout << "\n";
  cout << "  Reading the file " << file_in_name << ".\n";
  cout << flush;

  file_in.open ( file_in_name.c_str ( ) );

  if ( !file_in ) 
  {
    cout << "\n";
    cout << "TEST02 - Fatal error!\n";
    cout << "  Could not open the file " << file_in_name << "\n";
    return true;
  }

  cout << "\n";
  cout << "  The file " << file_in_name << " has been opened.\n";
  cout << flush;

  error = ppmb_read_header ( file_in, xsize, ysize, maxrgb );

  if ( error )
  {
    cout << "\n";
    cout << "  PPMB_READ_HEADER failed!\n";
    return error;
  }

  cout << "\n";
  cout << "  The header was read successfully.\n";
  cout << "  Number of rows of data    = " << xsize << "\n";
  cout << "  Number of columns of data = " << ysize << "\n";
  cout << "  Maximum data value =        " << maxrgb << "\n";
  cout << flush;

  r = new unsigned char[ xsize * ysize ];
  g = new unsigned char[ xsize * ysize ];
  b = new unsigned char[ xsize * ysize ];

  error = ppmb_read_data ( file_in, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "  PPMB_READ_DATA failed!\n";
    return error;
  }

  cout << "\n";
  cout << "  The data was read successfully.\n";

  delete [] r;
  delete [] g;
  delete [] b;

  return error;
}
//****************************************************************************80

void timestamp ( void )

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
