# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int file_line_count ( char input_filename[] );
void handle ( char input_filename[], int lines_per_page );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PAGES.
//
//  Discussion:
//
//    PAGES counts the lines in one or more files, and divides by 60.
//
//  Usage:
//
//      pages file(s)
//
//    or
//
//      pages
//
//    in which case the program will prompt the user for (one) filename.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2003
//
//  Author:
//
//    John Burkardt
//
{
  char input_filename[80];
  int i;
  int lines_per_page = 60;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "PAGES:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Estimate the number of pages in one or more files.\n";
    cout << "  A page is estimated to have " << lines_per_page << " lines.\n";
    cout << "\n";
  }
//
//  There must be an input file specified.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PAGES:\n";
    cout << "  Please enter an input file name:\n";

    cin.getline ( input_filename, sizeof ( input_filename ) );

    handle ( input_filename, lines_per_page );
  }
  else
  {
    for ( i = 1; i < argc; ++i )
    {

      strcpy ( input_filename, argv[i] );
      handle ( input_filename, lines_per_page );

    }
  }

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "PAGES:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void handle ( char input_filename[], int lines_per_page )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE processes one file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char INPUT_FILENAME[], the name of the file.
//
//    Input, int LINES_PER_PAGE, the number of lines per page.
//
{
  int lines;
  int pages;

  lines = file_line_count ( input_filename );

  pages = 1 + ( lines - 1 ) / lines_per_page;

  if ( pages <= 0 )
  {
    cout << "No pages in \"" << input_filename << "\".\n";
  } 
  else if ( pages == 1 )
  {
    cout << "1 page of " 
         << lines_per_page << " lines per page in \"" 
         << input_filename << "\".\n";
  } 
  else
  {
    cout << pages << " pages of " 
         << lines_per_page << " lines per page in \"" 
         << input_filename << "\".\n";
  }
  return;
}
//****************************************************************************80

int file_line_count ( char input_filename[] )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_LINE_COUNT counts the lines in one file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char INPUT_FILENAME[], the name of the file.
//
//    Output, int FILE_LINE_COUNT, the number of lines in the file.
//    This is set to -1 if the file cannot be opened.
//
{
  char c;
  ifstream input_unit;
  int lines;

  lines = 0;
//
//  Open the input file.
//
  input_unit.open ( input_filename );

  if ( !input_unit )
  {
    cout << "\n";
    cout << "FILE_LINE_COUNT:\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    return -1;
  }
//
//  Count line returns.
//
  while ( 1 )
  {
    input_unit.get ( c );

    if ( input_unit.eof() )
    {
      break;
    }

    if ( c == '\n' )
    {
      lines = lines + 1;
    }

  }
//
//  Close the file.
//
  input_unit.close ( );

  return lines;
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
