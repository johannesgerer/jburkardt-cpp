# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
void handle ( char input_file_name[], int *wide_line_width, 
int *wide_line_number );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WIDTH.
//
//  Discussion:
//
//    WIDTH reports the length of the longest line in each command line file.
//
//    This program reads every line of a file (terminated by a carriage return)
//    and reports the length of the longest one.
//
//  Usage:
//
//    width file(s)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 April 2003
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  char input_file_name[80];
  bool VERBOSE = false;
  int wide_line_number;
  int wide_line_width;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "WIDTH:\n";
    cout << "  C++ version\n";
    cout << "  Report length of longest line in file.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "WIDTH:\n";
    cout << "  Please enter the input file name:\n";

    cin.getline ( input_file_name, sizeof ( input_file_name ) );

    handle ( input_file_name, &wide_line_width, &wide_line_number );

    cout << "  The longest line of \"" << input_file_name 
         << "\" has length " << wide_line_width 
         << " and occurs at position " << wide_line_number << ".\n";
  }
//
//  Otherwise, get the file(s) from the argument list. 
//
  else 
  {

    for ( i = 1 ; i < argc ; ++i ) 
    {
      handle ( argv[i], &wide_line_width, &wide_line_number );

      cout << "  The longest line of \"" << argv[i] 
           << "\" has length " << wide_line_width 
           << " and occurs at position " << wide_line_number << ".\n";
    }

  } 

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "WIDTH:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void handle ( char input_file_name[], int *wide_line_width, 
  int *wide_line_number )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE reports the length of the longest line in a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char INPUT_FILE_NAME[], the name of the input file.
//
//    Output, int *WIDE_LINE_WIDTH, the length of the longest line.
//
//    Output, int *WIDE_LINE_NUMBER, the location of the longest line.
{
  int big_number;
  int big_width;
  char c;
  ifstream input_file;
  int input_file_width;
  int line_number;
  int line_width;
  bool VERBOSE = false;

  big_width = -1;
  big_number = -1;
//
//  Open the input file.
//
  input_file.open( input_file_name );

  if ( !input_file )
  {
    cout << "\n";
    cout << "WIDTH - Fatal error!\n";
    cout << "  Cannot open the input file " << input_file << ".\n";
    return;
  }
//
//  Examine characters.
//
  big_width = 0;
  line_width = 0;
  line_number = 0;

  while ( 1 ) 
  {
    input_file.get ( c );
 
    if ( input_file.eof ( ) )
    {
      break;
    }

    if ( c == '\n' ) 
    {
      line_number = line_number + 1;

      if ( big_width < line_width )
      {
        big_width = line_width;
        big_number = line_number;
      }
      line_width = 0;
    }
    else
    {
      line_width = line_width + 1;
    }

  }
//
//  Close the files.
//
  input_file.close ( );

  *wide_line_width = big_width;
  *wide_line_number = big_number;

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
