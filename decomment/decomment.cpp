# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void handle ( string input_filename, string output_filename );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DECOMMENT.
//
//  Discussion:
//
//    DECOMMENT makes a copy of a file, skipping all comment lines.
//
//    A comment line is exactly any line that begins with the "#" character.
//    I want this program because I commonly use a data file format that includes
//    such comment lines, but occasionally I need to delete such comments in
//    order to make an input file acceptable to a more uncouth program.
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
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "DECOMMENT:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Make a copy of a file, in which all # comments are gone.\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "DECOMMENT:\n";
    cout << "  Please enter the name of the file to be decommented.\n";
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
    cout << "DECOMMENT:\n";
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
    cout << "DECOMMENT:\n";
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
//    The routine reads a line, and if it does not begin with '#' it copies
//    it to the output file.
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
{
  int comment_count;
  ifstream input;
  ofstream output;
  string line;
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

  comment_count = 0;

  while ( 1 )
  {
    getline ( input, line );

    if ( line[0] != '#' )
    {
      output << line << "\n";
    }
    else
    {
      comment_count = comment_count + 1;
    }
    if ( input.eof() )
    {
      break;
    }
  }
//
// Close the files.
//
  input.close ( );
  output.close ( );

  cout << input_filename << " contained " << comment_count << " comments.\n";

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
