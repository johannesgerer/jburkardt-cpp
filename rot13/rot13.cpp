# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_to_rot13 ( char ch );
void handle ( string input_filename, string output_filename );
string s_to_rot13 ( string s );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ROT13.
//
//  Discussion:
//
//    ROT13 applies the "rot13" transformation to a file.
//
//    This program reads a file, and writes a copy in which characters have
//    been "rotated" 13 positions and digits have been "rotated" 5 positions.
//
//  Usage:
//
//    rot13 file1.txt
//
//    creates a rotated copy named "svyr6.gkg"
//
//    rot13 svyr5.gkg
//
//    recovers the original file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2011
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
    cout << "ROT13:\n";
    cout << "  C++ version\n";
	cout << "\n";
    cout << "  Copy a file, applying the rot13 rotation to characters,\n";
	cout << "  as well as a rot5 rotation to digits.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not on the command line, get it now.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "ROT13:\n";
    cout << "  Please enter the INPUT file name:\n";
    cin >> input_filename;
  }
  else 
  {
    input_filename = argv[1];
  }
//
//  Create the output file name.
//
  output_filename = s_to_rot13 ( input_filename );
  
  cout << "\n";
  cout << "  Output file name is \"" << output_filename << "\"\n";
//
//  Now we know the input and output file names, so go to it.
//
  handle ( input_filename, output_filename );
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "ROT13:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

char ch_to_rot13 ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_ROT13 converts a character to its ROT13 equivalent.
//
//  Discussion:
//
//    Two applications of CH_TO_ROT13 to a character will return the original.!
//
//    As a further scrambling, digits are similarly rotated using
//    a "ROT5" scheme.
//
//  Example:
//
//    Input:  Output:
//
//    a       n
//    C       P
//    J       W
//    1       6
//    5       0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, character CH, the character to be converted.
//
//    Output, character CH_TO_ROT13, the ROT13 equivalent of the character.
//
{
  char rot13;
//
//  [0:4] -> [5:9]
//
  if ( '0' <= ch && ch <= '4' )
  {
    rot13 = ch + 5;
  }
//
//  [5:9] -> [0:4]
//
  else if ( '5' <= ch && ch <= '9' )
  {
    rot13 = ch - 5;
  }
//
//  [A:M] -> [N:Z]
//
  else if ( 'A' <= ch && ch <= 'M' )
  {
    rot13 = ch + 13;
  }
//
//  [N:Z] -> [A:M]
//
  else if ( 'N' <= ch && ch <= 'Z' )
  {
    rot13 = ch - 13;
  }
//
//  [a:m] -> [n:z]
//
  else if ( 'a' <= ch && ch <= 'm' )
  {
    rot13 = ch + 13;
  }
//
//  [n:z] -> [a:m]
//
  else if ( 'n' <= ch && ch <= 'z' )
  {
    rot13 = ch - 13;
  }
  else
  {
    rot13 = ch;
  }

  return rot13;
}
//****************************************************************************80

void handle ( string input_filename, string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE makes a copy of a file after applying ROT13 to it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2011
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
  char c;
  int char_num;
  ifstream input;
  ofstream output;
  int line_num;
  bool VERBOSE = true;
//
//  Open the input and output files.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "ROT13::HANDLE - Fatal error!\n";
    cout << "  Cannot open the input file " << input_filename << ".\n";
    exit ( 1 );
  }

  output.open ( output_filename.c_str ( ) );

  if ( !output ) 
  {
    cout << "\n";
    cout << "ROT13::HANDLE - Fatal error!\n";
    cout << "  Cannot open the output file " << output_filename << ".\n";
    exit ( 1 );
  }
//
//  Transfer characters from the input file to the output file.
//
  char_num = 0;
  line_num = 0;

  while ( 1 ) 
  {
//
//  Get a new character.
//
    input.get ( c );

    if ( input.eof ( ) )
    {
      break;
    }

    char_num = char_num + 1;

    if ( c == '\n' ) 
    {
      line_num = line_num + 1;
    }
	
	c = ch_to_rot13 ( c );

    output.put ( c );
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
    cout << "    " <<  char_num  << " characters.\n";
    cout << "    " <<  line_num  << " lines.\n";
  }

  return;
}
//****************************************************************************80

string s_to_rot13 ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_ROT13 "rotates" the alphabetical characters in a string by 13 positions.
//
//  Discussion:
//
//    Two applications of the routine will return the original string.
//
//  Example:
//
//    Input:                      Output:
//
//    abcdefghijklmnopqrstuvwxyz  nopqrstuvwxyzabcdefghijklm
//    Cher                        Pure
//    James Thurston Howell       Wnzrf Guhefgba Ubjryy
//    0123456789                  5678901234
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//   29 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be "rotated".
//
//    Output, string S_TO_ROT13, the rotated string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    s2[i] = ch_to_rot13 ( s2[i] );
  }
  return s2;
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
