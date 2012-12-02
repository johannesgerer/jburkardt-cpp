# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int handle ( char file_in_name[], char file_out_name[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BITS_TO_UI.
//
//  Discussion:
//
//    BITS_TO_UI is a C++ program that converts a text file of
//   '0' and '1' characters to a binary file of 32 bit integers.
//
//    The binary file output by BITS_TO_UI is suitable for input
//    to Marsaglia's randomness testing program DIEHARD.
//
//  Usage:
//
//    bits_to_ui file1 file2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2003
//
//  Author:
//
//    John Burkardt
//
{
  char file_in_name[80];
  char file_out_name[80];
  int result;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "BITS_TO_UI:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Read an input text file of 0's and 1's.\n";
    cout << "  Write a binary output file in which groups of 32 bits are\n";
    cout << "  interpreted as unsigned integers.\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "BITS_TO_UI:\n";
    cout << "  Please enter the INPUT file name:\n";

    cin.getline ( file_in_name, sizeof ( file_in_name ) );
  }
  else
  {
    strcpy ( file_in_name, argv[1] );
  }
//
//  If the output file was not on the command line, get it now.
//
  if ( argc < 3 )
  {
    cout << "\n";
    cout << "BITS_TO_UI:\n";
    cout << "  Please enter the OUTPUT file name:\n";

    cin.getline ( file_out_name, sizeof ( file_out_name ) );
  }
  else
  {
    strcpy ( file_out_name, argv[2] );
  }
//
//  Now we know the input and output file names, so go to it.
//
  result = handle ( file_in_name, file_out_name );
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "BITS_TO_UI:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return result;
}
//****************************************************************************80

int handle ( char file_in_name[], char file_out_name[] )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE processes a given input file to an output file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char FILE_IN_NAME[], the name of the input file.
//
//    Input, char FILE_OUT_NAME[], the name of the output file.
//
//    Output, int HANDLE, is 0 for success, or 1 if there was an error.
//
{
  int bit_num;
  char c;
  ifstream file_in;
  ofstream file_out;
  int in_char;
  int out_ui;
  unsigned int ui;
  bool VERBOSE = true;
//
//  Open the input and output files.
//
  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open the input file " << file_in_name << ".\n";
    return 1;
  }

  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open the output file " << file_out_name << ".\n";
    return 1;
  }
//
//  Transfer characters from the input file to the output file.
//
  in_char = 0;
  out_ui = 0;
  bit_num = 0;
  ui = 0;

  while ( 1 )
  {
    file_in.get ( c );

    if ( file_in.eof ( ) )
    {
      break;
    }

    if ( c != '0' && c != '1' )
    {
      cout << "\n";
      cout << "HANDLE - Fatal error!\n";
      cout << "  Illegal input character C = " << c << "\n";
      cout << "  Input characters must be '0' or '1' only!\n";
      return 1;
    }

    in_char = in_char + 1;
    bit_num = bit_num + 1;

    ui = 2 * ui;

    if ( c == '1' )
    {
      ui = ui + 1;
    }
//
//  I don't understand why the simple act of writing a single unsigned integer
//  to a binary file requires the convoluted syntax below, but that's what
//  the manual says!
//
    if ( bit_num == 32 )
    {
      file_out.write ( ( char * )( &ui ), sizeof ( ui ) );
      out_ui = out_ui + 1;
      bit_num = 0;
      ui = 0;
    }
  }
//
//  Close the files.
//
  file_in.close ( );

  file_out.close ( );
//
//  Report.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "  The input file contains  " << in_char << " characters.\n";
    cout << "  The output file contains " << out_ui  << " unsigned integers.\n";
    cout << "  The number of leftover bits was " << bit_num << "\n";
  }

  return 0;
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
//    02 October 2003
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
