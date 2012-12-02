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
//    DETROFF removes obnoxious "character" + backspace pairs from a file. 
//
//  Usage:
//
//    detroff file1 file2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2009
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
    cout << "DETROFF:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Remove all character+backspace pairs in a file.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "DETROFF:\n";
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
    cout << "DETROFF:\n";
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
  handle ( input_filename, output_filename );

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "DETROFF:\n";
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
//    HANDLE makes a copy of a file minus "character"+backspace's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, string OUTPUT_FILENAME, the name of the output file.
//
{
  int bs_num;
  char c1;
  char c2;
  int cr_num;
  bool c2_saved;
  ifstream input;
  int input_num;
  int line_num;
  ofstream output;
  int output_num;
  bool VERBOSE = false;
//
//  Open the input and output files.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "DETROFF - Fatal error!\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    return;
  }

  output.open ( output_filename.c_str ( ) );

  if ( !output ) 
  {
    cout << "\n";
    cout << "DETROFF - Fatal error!\n";
    cout << "  Cannot open the output file " << output_filename << ".\n";
    return;
  }
//
//  Transfer characters from the input file to the output file.
//
  input_num = 0;
  output_num = 0;
  line_num = 0;
  bs_num = 0;
  c2_saved = false;

  while ( 1 ) 
  {
    input.get ( c1 );

    if ( input.eof ( ) )
    {
      if ( c2_saved )
      {
        output_num = output_num + 1;
        output.put ( c2 );
        c2_saved = false;
      }
      break;
    }

    input_num = input_num + 1;

    if ( c1 == '\n' ) 
    {
      if ( c2_saved )
      {
        output_num = output_num + 1;
        output.put ( c2 );
        c2_saved = false;
      }
      line_num = line_num + 1;
      output_num = output_num + 1;
      output.put ( c1 );
    }
    else if ( c1 == '\b' )
    {
      bs_num = bs_num + 1;
      c2_saved = false;
    }
    else if ( c2_saved )
    {
      output_num = output_num + 1;
      output.put ( c2 );
      c2 = c1;
    }
    else
    {
      c2_saved = true;
      c2 = c1;
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
    cout << "  The input file contains  " << input_num  << " characters.\n";
    cout << "  The input file contains  " << bs_num     << " backspaces.\n";
    cout << "  The output file contains " << output_num << " characters.\n";
    cout << "  Both files contain       " << line_num   << " lines.\n";
  }

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
