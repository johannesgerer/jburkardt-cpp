# include <cstdlib>
# include <iostream>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cctype>
# include <cstring>

using namespace std;

# include "filum.hpp"

//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  return ( ch1 == ch2 );
}
//****************************************************************************80

bool ch_is_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to be analyzed.
//
//    Output, bool CH_IS_DIGIT, is TRUE if the character is a digit.
//
{
  if ( '0' <= ch && ch <= '9' )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80*

char ch_low ( char ch )

//****************************************************************************80*
//
//  Purpose:
//
//    CH_LOW lowercases a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "tolower" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to lowercase.
//
//    Output, char CH_LOW, the lowercase character.
//
{
  if ( 65 <= ch && ch <= 90 ) 
  {
    ch = ch + 32;
  }   

  return ch;
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
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
//    Two applications of CH_TO_ROT13 to a character will return the original.
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

char digit_inc ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_INC increments a decimal digit.
//
//  Example:
//
//    Input  Output
//    -----  ------
//    '0'    '1'
//    '1'    '2'
//    ...
//    '8'    '9'
//    '9'    '0'
//    'A'    'A'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, a digit to be incremented.
//
//    Output, char DIGIT_INC, the incremented digit.
//
{
  if ( '0' <= c && c <= '8' )
  {
    return ( c + 1 );
  }
  else if ( c == '9' ) 
  {
    return '0';
  }
  else
  {
    return c;
  }
}
//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'  
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

int file_char_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_CHAR_COUNT counts the number of characters in a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME the name of the file.
//
//    Output, int FILE_CHAR_COUNT, the number of characters in the file.
//
{
  char c;
  ifstream input;
  int nchar;

  nchar = 0;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    nchar = -1;
    cerr << "\n";
    cerr << "FILE_CHAR_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return nchar;
  }
//
//  Read one character at a time.
//
  for ( ; ; )
  {
    input.read ( &c, 1 );
    
    if ( input.bad ( ) )
    {
      break;
    }
    if ( input.gcount ( ) != 1 )
    {
      break;
    }
    nchar = nchar + 1;
  }

  input.close ( );

  return nchar;
}
//****************************************************************************80

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words, 
//    separated by spaces.  There may also be some blank lines, and some 
//    comment lines, which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
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
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  string text;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) <= 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input >> text;

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
      {
        continue;
      }
      got_one = true;
      break;
    }
  }

  input.close ( );

  if ( !got_one )
  {
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

int file_delete ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_DELETE deletes a named file if it exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_DELETE, is 0 if the file deletion was successful.
//
{
  int value;
//
//  Does the file exist?
//
  if ( !file_exist ( filename ) )
  {
    return 1;
  }
//
//  Try to remove it.
//
  value = remove ( filename.c_str ( ) );

  if ( value != 0 )
  {
    cerr << "\n";
    cerr << "FILE_DELETE: Warning!\n";
    cerr << "  Could not delete \"" << filename << "\".\n";
    return value;
  }

  cout << "\n";
  cout << "FILE_DELETE:\n";
  cout << "  Deleting old version of \"" << filename << "\".\n";

  return value;
}
//****************************************************************************80

bool file_exist ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_EXIST reports whether a file exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;
  bool value;

  file.open ( filename.c_str ( ), ios::in );

  if ( !file )
  {
    value = false;
  }
  else
  {
    value = true;
  }
  return value;
}
//****************************************************************************80

string file_line_uniform ( string filename, int *seed, int *line_index, 
  int *line_num )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_LINE_UNIFORM returns a random line from a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    The algorithm used is interesting because it does not require
//    the number of lines in the file to be known in advance, and it
//    only reads the file once.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Tom Christiansen and Nathan Torkington,
//    "8.6: Picking a Random Line from a File",
//    Perl Cookbook, pages 284-285,
//    O'Reilly, 1999.
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int *LINE_INDEX, the index of the chosen line.
//
//    Output, int *LINE_NUM, the number of lines in the file.
//
//    Output, string FILE_LINE_UNIFORM, a random line from the file.
//
{
  ifstream input;
  double r;
  string text;
  string text2;

  *line_num = 0;
  *line_index = -1;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_LINE_UNIFORM - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    text2 = "";
    return text2;
  }
//
//  Read the lines.
//
  for ( ; ; )
  {
    getline ( input, text, '\n' );

    if ( input.eof ( ) )
    {
      break;
    }

    *line_num = *line_num + 1;

    r = r8_uniform_01 ( seed );

    if ( r * ( double ) ( *line_num ) <= 1.0 )
    {
      text2 = text;
      *line_index = *line_num;
    }

  }
  input.close ( );

  return text2;
}
//****************************************************************************80

int file_line_width ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_LINE_WIDTH reports the length of the longest line in a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_LINE_WIDTH, the length of the longest line.
//
{
  char c;
  ifstream input;
  int line_width;
  int value;

  value = -1;
//
//  Open the input file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_LINE_WIDTH - Fatal error!\n";
    cerr << "  Cannot open the input file \"" << filename << "\".\n";
    return value;
  }
//
//  Examine characters.
//
  value = 0;
  line_width = 0;

  while ( 1 ) 
  {
    input.get ( c );
 
    if ( input.eof ( ) )
    {
      break;
    }

    if ( c == '\n' || c == '\r' ) 
    {
      line_width = 0;
    }
    else
    {
      line_width = line_width + 1;
      value = i4_max ( value, line_width );
    }
  }
//
//  Close the file.
//
  input.close ( );

  return value;
}
//****************************************************************************80

int file_para_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_PARA_COUNT counts the number of paragraphs in a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.  A paragraph is
//    a sequence of nonblank lines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_PARA_COUNT, the number of paragraphs found in the file.
//
{
  ifstream input;
  int lenc;
  int lenc_old;
  int para_num;
  string text;

  para_num = 0;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    para_num = -1;
    cerr << "\n";
    cerr << "FILE_PARA_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    exit ( 1 );
  }
//
//  Read one line at a time.
//
  lenc = 0;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    lenc_old = lenc;
    lenc = s_len_trim ( text );

    if ( 0 < lenc && lenc_old <= 0 )
    {
      para_num = para_num + 1;
    }
  }

  input.close ( );

  return para_num;
}
//****************************************************************************80

int file_row_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int record_num;
  int row_num;
  string text;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the file: \"" << filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( text[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( text ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;
  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

int file_sequence_delete ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_SEQUENCE_DELETE deletes a file sequence.
//
//  Discussion:
//
//    We suppose the user has a set of files whose names differ only
//    in some numeric tag that is sequentially increasing, as, perhaps,
//    "file001.txt", "file002.txt" through "file137.txt", say.
//
//    The user specifies FILENAME as the name of the first file in the
//    sequence.  This function deletes that file, generates the next
//    name in the sequence, and, if a file with that name exists, it
//    deletes it as well.  The process continues until a file name is
//    reached for which there is no existing file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the first file
//    in the sequence.
//
//    Output, int FILE_SEQUENCE_DELETE, the number of files deleted.
//
{
  int delete_count;
  string filename2;

  delete_count = 0;

  filename2 = filename;

  while ( file_exist ( filename2 ) )
  {
    file_delete ( filename2 );

    delete_count = delete_count + 1;

    filename_inc ( &filename2 );
  }

  return delete_count;
}
//****************************************************************************80

void file_sequence_size ( string filename, int *file_dim, int *file_num )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_SEQUENCE_SIZE sizes a file sequence.
//
//  Discussion:
//
//    We suppose the user has a set of files whose names differ only
//    in some numeric tag that is sequentially increasing, as, perhaps,
//    "file001.txt", "file002.txt" through "file137.txt", say.
//
//    The user specifies the name of the first file in the sequence.
//    This function determines the number of files in the sequence, 
//    and makes a guess for the "dimension" of the files, that is, the number
//    of numeric data items.
//
//    Note that the function only checks the dimension of the data in
//    the first file.  It is up to the user to determine whether this
//    dimension is used for every file in the sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the first file in 
//    the sequence.
//
//    Output, int *FILE_DIM, the dimension of the data in one file.
//
//    Output, int *FILE_NUM, the number of files.
//
{
  string filename2;

  *file_num = 0;
  *file_dim = 0;

  filename2 = filename;

  for ( ; ; )
  {
    if ( !file_exist ( filename2 ) )
    {
      break;
    }
    *file_num = *file_num + 1;

    if ( *file_num == 1 )
    {
      *file_dim = file_word_count ( filename2 );
    }
    filename_inc ( &filename2 );
  }

  return;
}
//****************************************************************************80

int file_word_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_WORD_COUNT counts the number of words in a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_WORD_COUNT, the number of words found.
//
{
  char c;
  ifstream input;
  bool white;
  int word_num;

  word_num = 0;
  white = true;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    word_num = -1;
    cerr << "\n";
    cerr << "FILE_WORD_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return word_num;
  }
//
//  Read one character at a time.
//
  for ( ; ; )
  {
    input.read ( &c, 1 );
    
    if ( input.bad ( ) )
    {
      break;
    }
    if ( input.gcount ( ) != 1 )
    {
      break;
    }
//
//  A new word occurs every time we were reading whitespace but encounter 
//  a nonwhitespace character.
//
    if ( !isspace ( c ) && white )
    {
      word_num = word_num + 1;
    }

    white = isspace ( c );
  }

  input.close ( );

  return word_num;
}
//****************************************************************************80

void filename_dec ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_DEC decrements a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be decreased by 1 on
//    each call.  If this number is all 0's on input, the output number
//    is all 9's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty or contains no digits, an error condition occurs.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to12.txt"     "a7to11.txt"  (typical case.  Last digit decremented)
//      "a8to00.txt"     "a7to99.txt"  (last digit decremented, with carry.)
//      "a0to00.txt"     "a9to99.txt"  (wrap around)
//      "cat.txt"        " "           (no digits in input name)
//      " "              Stop!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the string to be decremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_DEC - Fatal error!\n";
    cerr << "  Input file name is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '0' )
      {
        c = '9';
        (*filename)[i] = c;
      }
      else
      {
        c = c - 1;
        (*filename)[i] = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_DEC - Fatal error!\n";
    cerr << "  Input file name contained no digits.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void filename_ext_get ( string filename, int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_EXT_GET determines the "extension" of a file name.
//
//  Discussion:
//
//    The "extension" of a file name is the string of characters
//    that appears after the LAST period in the name.  A file
//    with no period, or with a period as the last character
//    in the name, has a "null" extension.
//
//    Blanks are unusual in file names.  This routine ignores all
//    trailing blanks, but will treat initial or internal blanks
//    as regular characters acceptable in a file name.
//
//    If no period occurs in FILENAME, then
//      I = J = -1;
//    Otherwise,
//      I is the position of the LAST period in FILENAME, and J is the
//      position of the last nonblank character following the period.
//
//  Example:
//
//    FILENAME    I  J
//
//    bob.for     3  6
//    N.B.C.D     5  6
//    Naomi.      5  5
//    Arthur     -1 -1
//    .com        0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, a file name to be examined.
//
//    Output, int *I, *J, the indices of the first and last characters
//    in the file extension.
//
{
  size_t location;

  location = filename.rfind ( "." );

  if ( location == string::npos )
  {
    *j = -1;
  }
  else
  {
    *j = filename.length ( ) - 1;
  }
  *i = ( int ) location;

  return;
}
//****************************************************************************80

string filename_ext_swap ( string filename, string ext )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_EXT_SWAP replaces the current "extension" of a file name.
//
//  Discussion:
//
//    The "extension" of a file name is the string of characters
//    that appears after the LAST period in the name.  A file
//    with no period, or with a period as the last character
//    in the name, has a "null" extension.
//
//  Example:
//
//          Input           Output
//    ================     ==================
//    FILENAME     EXT     FILENAME_EXT_SWAP
//
//    bob.for      obj     bob.obj
//    bob.bob.bob  txt     bob.bob.txt
//    bob          yak     bob.yak
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, a file name.
//
//    Input, string EXT, the extension to be added to the file name.
//
//    Output, string FILENAME_EXT_SWAP, the file name with the new extension.
//
{
  string filename2;
  size_t i;
//
//  Look for the LAST occurrence of a period.
//
  i = filename.rfind ( "." );

  if ( i == string::npos ) 
  {
    filename2 = filename + "." + ext;

  }
  else
  {
    filename2 = filename.substr ( 0, i + 1 ) + ext;
  }

  return filename2;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

void filename_inc_nowrap ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC_NOWRAP increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the (nonempty) name contains no digits, or all the digits are
//    9, then the empty string is returned.
//
//    If the empty string is input, the routine stops.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a8to99.txt"     "a9to00.txt"
//      "a9to99.txt"     " "
//      "cat.txt"        " "
//      " "              STOP!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int carry;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC_NOWRAP - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;
  carry = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      carry = 0;

      if ( c == '9' )
      {
        carry = 1;
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  Unsatisfied carry.  The input digits were all 9.  Return blank.
//
  if ( carry == 1 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1 and I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  if ( i2 < i1 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1 and I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  if ( i1 < i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 ) 
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

void number_inc ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    NUMBER_INC increments the integer represented by a string.
//
//  Discussion:
//
//    If the string contains characters that are not digits, they will
//    simply be ignored.  If the integer is all 9's on input, then
//    the output will be all 0's.
//
//  Example:
//
//    Input      Output
//    -----      ------
//    '17'       '18'
//    'cat3'     'cat4'
//    '2for9'    '3for0'
//    '99thump'  '00thump'
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
//    Input/output, string S, a string representing an integer.
//
{
  int i;
  int s_length;

  s_length = s.length ( );

  for ( i = s_length - 1; 0 <= i; i-- )
  {
    if ( ch_is_digit ( s[i] ) )
    {
      s[i] = digit_inc ( s[i] );

      if ( s[i] != '0' )
      {
        return;
      }
    }
  }

  return;
}
//*********************************************************************

float r4_abs ( float x )

//*********************************************************************
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         R8_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
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
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;
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

void s_low ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LOW lowercases a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be lowercased.  
//    On output, the string has been lowercased.
//
{
  size_t i;

  for ( i = 0; i < s.length ( ); i++ ) 
  {
    s[i] = ch_low ( s[i] );
  }

  return;
}
//****************************************************************************80

string s_reverse ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_REVERSE reverses the characters in a string.
//
//  Example:
//
//    Input        Output
//
//    ' Cat'       'taC '
//    'Goo gol  '  'log ooG  '
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to reverse.
//
//    Output, string S_REVERSE, the reversed string.
//
{
  int i;
  int ihi;
  int j;
  int s_length;
  string s2;
  char t;

  s_length = s.length ( );
  s2 = s;

  ihi = ( s_length / 2 ) - 1;

  for ( i = 0 ; i <= ihi; i++ )
  {
    j = s_length - i - 1;
    t     = s2[i];
    s2[i] = s2[j];
    s2[j] = t;
  }

  return s2;
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

bool s_to_i4vec ( string s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
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
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
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
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
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
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

void s_to_rot13 ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_ROT13 "rotates" the characters in a string by 13 positions.
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
//   01 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string S, a string to be "rotated".
//
{
  int i;
  int n;

  n = s.length ( );

  for ( i = 0; i < n; i++ )
  {
    s[i] = ch_to_rot13 ( s[i] );
  }

  return;
}
//****************************************************************************80

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
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
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

