# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

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
//    MAIN is the main program for UNCONTROL.
//
//  Discussion:
//
//    UNCONTROL makes a copy of a file with all control characters removed.
//
//    The one exception is the NEWLINE character, which is not removed.
//
//  Usage:
//
//    uncontrol
//     in which case the input and output files will be prompted for;
//
//    uncontrol file1
//      in which case the output file will be prompted for;
//
//    uncontrol file1 file2
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
    cout << "UNCONTROL:\n";
    cout << "  C++ version\n";
    cout << "  Remove all control characters (except newlines).\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "UNCONTROL:\n";
    cout << "  Please enter the name of the file to be filtered.\n";
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
    cout << "UNCONTROL:\n";
    cout << "  Please enter the OUTPUT file name:\n";
    cin >> output_filename;
  }
  else 
  {
    output_filename = argv[2];
  }

  result = handle ( input_filename, output_filename );
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "UNCONTROL:\n";
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
//    HANDLE makes a copy of a file in which CR's are replace by LF's.
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
  int control_num;
  ifstream input;
  ofstream output;
  int input_char;
  int output_char;
  bool VERBOSE = false;
//
//  Open the input and output files.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "UNCONTROL - Fatal error!\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    return 1;
  }

  output.open ( output_filename.c_str ( ) );

  if ( !output ) 
  {
    cout << "\n";
    cout << "UNCONTROL - Fatal error!\n";
    cout << "  Cannot open the output file " << output_filename << ".\n";
    return 1;
  }
//
//  Transfer characters from the input file to the output file.
//
  input_char = 0;
  output_char = 0;
  control_num = 0;

  while ( 1 ) 
  {
    input.get ( c );

    if ( input.eof ( ) )
    {
      break;
    }

    input_char = input_char + 1;

    if ( c == '\n' || ! iscntrl ( c ) )
    {
      output.put ( c );
      output_char = output_char + 1;
    }
    else
    {
      control_num = control_num + 1;
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
    cout << "  The input file contains  " <<  input_char << " characters.\n";
    cout << "  The output file contains " << output_char << " characters.\n";
    cout << "  The number of deleted control characters was " 
      << control_num << ".\n";
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
//    04 October 2003
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
