# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
int handle ( string input_filename, string output_filename );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    CR2CRLF replaces carriage returns by carriage returns + line feeds.
//
//  Usage:
//
//    cr2crlf
//      in which case the input and output files will be prompted for;
//
//    cr2crlf file1
//      in which case the output file will be prompted for;
//
//    cr2crlf file1 file2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  string input_filename;
  string output_filename;
  int result;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "CR2CRLF:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Replace carriage returns by carriage returns + line feeds.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "CR2CRLF:\n";
    cout << "  Please enter the INPUT file name:\n";
    cin >> input_filename;
  }
  else 
  {
    input_filename = argv[1];
  }
//
//  If the output file was not on the command line, get it now.
//
  if ( argc < 3 ) 
  {
    cout << "\n";
    cout << "CR2CRLF:\n";
    cout << "  Please enter the OUTPUT file name:\n";
    cin >> output_filename;
  }
  else 
  {
    output_filename = argv[2];
  }
//
//  Now we know the input and output file names, so go to it.
//
  result = handle ( input_filename, output_filename );

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "CR2CRLF:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return result;
}
//****************************************************************************80

int handle ( string input_filename, string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE makes a copy of a file in which CR's are replaced by CR+LF's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILE, the name of the input file.
//
//    Input, string OUTPUT_FILE, the name of the output file.
//
//    Output, int HANDLE, is 0 for success, or 1 if there was an error.
//
{
  char c;
  int cr_num;
  ifstream input;
  int input_char;
  ofstream output;
  int output_char;
  bool VERBOSE = false;
//
//  Open the input and output files.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "CR2CRLF - Fatal error!\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    return 1;
  }

  output.open ( output_filename.c_str ( ) );

  if ( !output ) 
  {
    cout << "\n";
    cout << "CR2CRLF - Fatal error!\n";
    cout << "  Cannot open the output file " << output_filename << ".\n";
    return 1;
  }
//
//  Transfer characters from the input file to the output file.
//
  input_char = 0;
  output_char = 0;
  cr_num = 0;

  while ( 1 ) 
  {
    input.get ( c );

    if ( input.eof ( ) )
    {
      break;
    }

    input_char = input_char + 1;

    if ( c == '\r' ) 
    {
      cr_num = cr_num + 1;
      output.put ( '\r' );
      output.put ( '\n' );
      output_char = output_char + 2;
    }
    else
    {
      output.put ( c );
      output_char = output_char + 1;
    }

  }
//
//  Close the files.
//
  input.close ( );

  output.close ( );
//
//  Report.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "  The input file contains  " << input_char  << " characters.\n";
    cout << "  The input file contains  " << cr_num      << " carriage returns.\n";
    cout << "  The output file contains " << output_char << " characters.\n";
  }

  return 0;
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
//    03 October 2003
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
