# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int handle ( string input_filename, string output_filename, int word_max );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WRAP2.
//
//  Discussion:
//
//    WRAP2 rewrites a file to have a given number of "words" per line.
//
//    This program reads a file, and writes a copy for which each line
//    has a fixed number of "words", that is, blank-delimited strings.
//
//  Usage:
//
//    wrap2 file1 file2 char_max
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
  bool error;
  string input_filename;
  int last;
  string output_filename;
  int result;
  string s;
  bool VERBOSE = false;
  int char_max = 80;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "WRAP2:\n";
    cout << "  C++ version\n";
    cout << "  Insert new lines where necessary into a copy of the input\n";
    cout << "  file, so that the lines of the output file are no longer\n";
    cout << "  than a given number of characters.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "WRAP2:\n";
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
    cout << "WRAP2:\n";
    cout << "  Please enter the OUTPUT file name:\n";
    cin >> output_filename;
  }
  else
  {
    output_filename = argv[2];
  }
//
//  If the number of words per line was not on the command line, get it now.
//
  if ( argc < 4 )
  {
    cout << "\n";
    cout << "WRAP2:\n";
    cout << "  Please enter the maximum number of characters per line:\n";
    cin >> char_max;
  }
  else
  {
    s = argv[3];
    char_max = s_to_i4 ( s, &last, &error );
  }
//
//  Now we know the input and output file names, so go to it.
//
  result = handle ( input_filename, output_filename, char_max );
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "WRAP2:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return result;
}
//****************************************************************************80

int handle ( string input_filename, string output_filename, int char_max )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE makes a copy of a file with a given number of characters per line.
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
//    Input, int CHAR_MAX, the maximum number of characters per line.
//
//    Output, int HANDLE, is 0 for success, or 1 if there was an error.
//
{
  int break_early;
  char c;
  char c_old;
  ifstream input;
  int input_char_num;
  int input_line_num;
  int line_char_num;
  ofstream output;
  int output_char_num;
  int output_line_num;
  bool VERBOSE = true;
//
//  Open the input and output files.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "WRAP2::HANDLE - Fatal error!\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    return 1;
  }

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "WRAP2::HANDLE - Fatal error!\n";
    cout << "  Cannot open the output file " << output_filename << ".\n";
    return 1;
  }
//
//  Transfer characters from the input file to the output file.
//
  input_char_num = 0;
  input_line_num = 0;
  line_char_num = 0;
  output_char_num = 0;
  output_line_num = 0;
  c = ' ';

  if ( 25 < char_max )
  {
    break_early = char_max - 12;
  }
  else
  {
    break_early = char_max / 2;
  }

  while ( 1 )
  {
//
//  Get a new character.
//
    c_old = c;

    input.get ( c );

    if ( input.eof ( ) )
    {
      output.put ( '\n' );
      output_char_num = output_char_num + 1;
      output_line_num = output_line_num + 1;
      break;
    }

    input_char_num = input_char_num + 1;

    if ( c == '\n' )
    {
      input_line_num = input_line_num + 1;
      line_char_num = 0;
    }
//
//  If we can break, should we?
//
    if ( c == ' ' && break_early <= line_char_num )
    {
      output.put ( '\n' );
      output_char_num = output_char_num + 1;
      output_line_num = output_line_num + 1;
      line_char_num = 0;
    }
//
//  If we must break, do so.
//
    else if ( char_max - 1 <= line_char_num )
    {
      output.put ( c );
      output.put ( '-' );
      output.put ( '\n' );
      output_char_num = output_char_num + 3;
      output_line_num = output_line_num + 1;
      line_char_num = 0;
    }
//
//  Output the character.
//
    else
    {
      output.put ( c );
      output_char_num = output_char_num + 1;
      line_char_num = line_char_num + 1;
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
    cout << "  The input file \"" << input_filename << "\" contains:\n";
    cout << "    " <<  input_char_num  << " characters.\n";
    cout << "    " <<  input_line_num  << " lines.\n";
    cout << "\n";
    cout << "  The output file \"" << output_filename << "\" contains:\n";
    cout << "    " <<  output_char_num << " characters.\n";
    cout << "    " <<  output_line_num << " lines.\n";
  }

  return 0;
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

int s_to_i4 ( string s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
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
//    Input, string S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for ( ; ; )
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
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
