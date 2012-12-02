# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void handle ( string input_filename, string output_filename );
int s_len_trim ( string s );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DEBLANK.
//
//  Discussion:
//
//    DEBLANK makes a copy of a file, skipping all blank lines.
//
//    A blank line is any line that contains no characters except blanks.
//
//    If the final line of the input file isn't blank, but does not 
//    terminate with a NEWLINE character, the corresponding line in
//    the output file will have a terminating NEWLINE appended to it.
//    A previous version of this program would instead inadvertently
//    omit this final line from the output file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  string input_filename;
  string output_filename;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "DEBLANK:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Make a copy of a file, in which all blank lines are gone.\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "DEBLANK:\n";
    cout << "  Please enter the name of the file to be deblanked.\n";
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
    cout << "DEBLANK:\n";
    cout << "  Please enter the OUTPUT file name:\n";
    cin >> output_filename;
  }
  else 
  {
    output_filename = argv[2];
  }

  handle ( input_filename, output_filename );

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "DEBLANK:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void handle ( string input_filename, string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE handles a single file.
//
//  Discussion:
//
//    The routine reads a line, and if it is not a blank line, copies
//    it to the output file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2010
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
{
  int blank_count;
  ifstream input;
  string line;
  ofstream output;
//
// Open the input file.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input ) 
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open \"" << input_filename << "\"\n";
    return;
  }
//
// Open the output file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output ) 
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open \"" << output_filename << "\"\n";
    return;
  }

  blank_count = 0;

  while ( 1 )
  {
    getline ( input, line );

    if ( 0 < s_len_trim ( line ) )
    {
      output << line << "\n";
    }
    else
    {
      blank_count = blank_count + 1;
    }
    if ( input.eof ( ) )
    {
      break;
    }
  }
//
//  Close the files.
//
  input.close ( );
  output.close ( );

  cout << input_filename << " contains " << blank_count << " blank lines.\n";

  return;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
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
