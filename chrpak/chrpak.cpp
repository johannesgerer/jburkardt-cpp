# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <sstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "chrpak.hpp"

//****************************************************************************80

int a_to_i4 ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    A_TO_I4 returns the index of an alphabetic character.
//
//  Example:
//
//    CH  A_TO_I4
//
//    'A'   1
//    'B'   2
//    ...
//    'Z'  26
//    'a'  27
//    'b'  28
//    ...
//    'z'  52
//    '$'   0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, a character.
//
//    Output, int A_TO_I4, is the alphabetic index of the character,
//    between 1 and 26 if the character is a capital letter,
//    between 27 and 52 if it is lower case, and 0 otherwise.
//
{
  int value;

  if ( 'A' <= ch && ch <= 'Z' )
  {
    value = ( int ) ( ch - 'A' + 1 );
  }
  else if ( 'a' <= ch && ch <= 'z' )
  {
    value = ( int ) ( ch - 'a' + 26 + 1 );
  }
  else
  {
    value = 0;
  }
  return value;
}
//****************************************************************************80

int base_to_i4 ( char *s, int base )

//****************************************************************************80
//
//  Purpose:
//
//    BASE_TO_I4 returns the value of an integer represented in some base.
//
//  Discussion:
//
//    BASE = 1 is allowed, in which case we allow the digits '1' and '0',
//    and we simply count the '1' digits for the result.
//
//    Negative bases between -16 and -2 are allowed.
//
//    The base -1 is allowed, and essentially does a parity check on
//    a string of 1's.
//
//  Example:
//
//        Input      Output
//    -------------  ------
//         S   BASE       I
//    ------  -----  ------
//      '101'     2       5
//    '-1000'     3     -27
//      '100'     4      16
//   '111111'     2      63
//   '111111'    -2      21
//   '111111'     1       6
//   '111111'    -1       0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string.  The elements of S are
//    blanks, a plus or minus sign, and digits.  Normally, the digits
//    are representations of integers between 0 and |BASE-1|.  In the
//    special case of base 1 or base -1, we allow both 0 and 1 as digits.
//
//    Input, int BASE, the base in which the representation is given.
//    Normally, 2 <= BASE <= 16.  However, there are two exceptions.
//
//    Output, int BASE_TO_I4, the integer.
//
{
  char c;
  int i;
  int ichr;
  int idig;
  int isgn;
  int istate;
  int nchar;

  nchar = charstar_len_trim ( s );

  if ( base == 0 )
  {
    cerr << "\n";
    cerr << "BASE_TO_I4 - Serious error!\n";
    cerr << "  The input base is zero.\n";
    i = -1;
    return i;
  }

  if ( 16 < abs ( base ) )
  {
    cerr << "\n";
    cerr << "BASE_TO_I4 - Serious error!\n";
    cerr << "  The input base is greater than 16!\n";
    i = -1;
    return i;
  }

  i = 0;
  istate = 0;
  isgn = 1 ;
  ichr = 1;

  while ( ichr <= nchar )
  {
    c = s[ichr-1];
//
//  Blank.
//
    if ( c == ' ' )
    {
      if ( istate == 2 )
      {
        break;
      }
    }
//
//  Sign, + or -.
//
    else if ( c == '-' )
    {
      if ( istate != 0 )
      {
        break;
      }
      istate = 1;
      isgn = -1;
    }
    else if ( c == '+' )
    {
      if ( istate != 0 )
      {
        break;
      }
      istate = 1;
    }
    else
//
//  Digit?
//
    {
      idig = hex_digit_to_i4 ( c );

      if ( abs ( base ) == 1 && ( idig == 0 || idig == 1 ) )
      {
        i = base * i + idig;
        istate = 2;
      }
      else if ( 0 <= idig && idig < abs ( base ) )
      {
        i = base * i + idig;
        istate = 2;
      }
      else
      {
        cerr << "\n";
        cerr << "BASE_TO_I4 - Serious error!\n";
        cerr << "  Illegal digit = \"" << c << "\"\n";
        cerr << "  Conversion halted prematurely!\n";
        i = -1;
        return i;
      }
    }
    ichr = ichr + 1;
  }
//
//  Once we're done reading information, we expect to be in state 2.
//
  if ( istate != 2 )
  {
    cerr << "\n";
    cerr << "BASE_TO_I4 - Serious error!\n";
    cerr << "  Unable to decipher input!\n";
    i = -1;
    return i;
  }
//
//  Account for the sign.
//
  i = isgn * i;

  return i;
}
//****************************************************************************80

int binary_to_i4 ( string s )

//****************************************************************************80
/*
  Purpose:

    BINARY_TO_I4 converts a binary representation into an I4.

  Example:

        S        I

      '101'      5
    '-1000'     -8
        '1'      1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2013

  Author:

    John Burkardt

  Parameters:

    Input, string S, the binary representation.

    Output, int BINARY_TO_I4, the I4 whose representation was input.
*/
{
  char c;
  int i;
  int ichr;
  int isgn;
  int s_len;
  int state;

  s_len = s.length ( );

  i = 0;
  ichr = 0;
  state = 0;
  isgn = 1;

  while ( ichr < s_len )
  {
    c = s[ichr];
//
//  Blank.
//
    if ( c == ' ' )
    {
      if ( state == 2 )
      {
        state = 3;
      }
    }
//
//  Sign, + or -.
//
    else if ( c == '-' )
    {
      if ( state == 0 )
      {
        state = 1;
        isgn = -1;
      }
      else
      {
        state = -1;
      }
    }
    else if ( c == '+' )
    {
      if ( state == 0 )
      {
        state = 1;
      }
      else
      {
        state = -1;
      }
    }
//
//  Digit, 0 or 1.
//
    else if ( c == '1' )
    {
      i = 2 * i;
      i = i + 1;
      state = 2;
    }
    else if ( c == '0' )
    {
      i = 2 * i;
      state = 2;
    }
//
//  Illegal or unknown sign.
//
    else
    {
      cout << "\n";
      cout << "BINARY_TO_I4 - Serious error!\n";
      cout << "  Illegal digit = '" << c << "'.\n";
      cout << "  Conversion halted prematurely!\n";
      exit ( 1 );
    }

    if ( state == -1 )
    {
      cout << "\n";
      cout << "BINARY_TO_I4 - Serious error!\n";
      cout << "  Unable to decipher input!\n";
      exit ( 1 );
    }

    if ( 3 <= state )
    {
      break;
    }

    ichr = ichr + 1;
  }
//
//  Apply the sign.
//
  i = isgn * i;

  return i;
}
//****************************************************************************80

void byte_to_int ( unsigned char *bvec, unsigned int *ival )

//****************************************************************************80
//
//  Purpose:
//
//    BYTE_TO_INT converts 4 bytes into an unsigned integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned char *BVEC, is a pointer to a character string.
//    The contents of BVEC through BVEC+3 are the bytes of IVAL,
//    from high order to low.
//
//    Output, unsigned int IVAL, the integer represented by the bytes.
//
{
  int i;

  *ival = 0;

  for ( i = 0; i < 4; i++ )
  {
    *ival = *ival << 8;
    *ival = *ival + *bvec;
    bvec = bvec + 1;
  }
  return;
}
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

void ch_count_cvec_add ( int n, unsigned char cvec[], int count[256] )

//****************************************************************************80
//
//  Purpose:
//
//    CH_COUNT_CVEC_ADD adds a character vector to a character count.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, unsigned char CVEC[n], a vector of characters.
//
//    Input/output, int COUNT[256], the character counts.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    count[cvec[i]] = count[cvec[i]] + 1;
  }

  return;
}
//****************************************************************************80

void ch_count_file_add ( string file_name, int count[256] )

//****************************************************************************80
//
//  Purpose:
//
//    CH_COUNT_FILE_ADD adds characters in a file to a character count.
//
//  Discussion:
//
//    Note that, although C is "really" an unsigned char, it must be
//    declared an INT so that it can have the value EOF when FGETC
//    reaches the end of file.  Otherwise...it never does, and you'll
//    be sorry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to examine.
//
//    Output, int COUNT[256], the character counts.
//
{
  char c;
  int i;
  ifstream input;
//
//  Open the file.
//
  input.open ( file_name.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "CH_COUNT_FILE_ADD - Fatal error!\n";
    cerr << "  Cannot open the input file " << file_name << ".\n";
    exit ( 1 );
  }

  while ( 1 )
  {
    input.get ( c );

    if ( input.eof ( ) )
    {
      break;
    }
    i = ( int ) c;
    count[i] = count[i] + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void ch_count_init ( int count[256] )

//****************************************************************************80
//
//  Purpose:
//
//    CH_COUNT_INIT initializes a character count.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int COUNT[256], the character counts.
//
{
  int i;

  for ( i = 0; i <= 255; i++ )
  {
    count[i] = 0;
  }

  return;
}
//****************************************************************************80

void ch_count_print ( int count[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    CH_COUNT_PRINT prints a set of character counts.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int COUNT[256], the character counts.
//
//    Input, char *TITLE, a title to be printed.
//
{
  int i;
  double percent;
  int total;

  total = 0;
  for ( i = 0; i <= 255; i++ )
  {
    total = total + count[i];
  }

  cout << "\n";
  cout << "  " << title << "\n";
  cout << "\n";
  cout << "  Char  Percent  Count\n";
  cout << "\n";

  for ( i = 0; i <= 255; i++ )
  {
    if ( 0 < count[i] )
    {
      if ( total == 0 )
      {
        percent = 0.0;
      }
      else
      {
        percent = ( double ) ( 100 * count[i] ) / ( double ) ( total );
      }
      if ( 32 <= i && i <= 126 )
      {
        cout << "     "   << (char) i << "  "
             << setw(3) << percent << "  "
             << count[i] << "\n";
      }
      else
      {
        cout << "  #"     << setw(3) << i << "  "
             << setw(6) << percent << "  "
             << count[i] << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void ch_count_s_add ( unsigned char *s, int count[256] )

//****************************************************************************80
//
//  Purpose:
//
//    CH_COUNT_S_ADD adds a character string to a character histogram.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned char *S, a string to be examined.
//
//    Input/output, int COUNT[256], the character counts.
//
{
  while ( *s )
  {
    count[*s] = count[*s] + 1;
    *s++;
  }

  return;
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
//    13 September 2009
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
  bool value;

  if ( 97 <= ch1 && ch1 <= 122 )
  {
    ch1 = ch1 - 32;
  }
  if ( 97 <= ch2 && ch2 <= 122 )
  {
    ch2 = ch2 - 32;
  }

  value = ( ch1 == ch2 );

  return value;
}
//****************************************************************************80

int ch_index_first ( string s, char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_INDEX_FIRST finds the first occurrence of a character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be searched.
//
//    Input, char C, the character to be searched for.
//
//    Output, int CH_INDEX_FIRST, the index of the first occurrence
//    of the character, or -1 if it does not occur.
//
{
  int i;
  int nchar;

  nchar = s.length ( );

  for ( i = 0; i < nchar; i++ )
  {
    if ( s[i] == c )
    {
      return i;
    }
  }

  return -1;
}
//****************************************************************************80

int ch_index_last ( string s, char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_INDEX_LAST finds the last occurrence of a character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be searched.
//
//    Input, char C, the character to be searched for in s.
//
//    Output, int CH_INDEX_LAST, the index of the last occurrence
//    of the character, or -1 if it does not occur.
//
{
  int i;
  int j;
  int nchar;

  j = - 1;

  nchar = s.length ( );

  for ( i = 0; i < nchar; i++ )
  {
    if ( s[i] == c )
    {
      j = i;
    }
  }

  return j;
}
//****************************************************************************80

bool ch_is_alpha ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_ALPHA is TRUE if a charaacter is alphabetic.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, a character to check.
//
//    Output, bool CH_IS_ALPHA is TRUE if the character is alphabetic.
//
{
  bool value;

  if ( ( 'a' <= ch && ch <= 'z' ) ||
       ( 'A' <= ch && ch <= 'Z' ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

bool ch_is_alphanumeric ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_ALPHANUMERIC is TRUE if a character is alphanumeric.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, a character to check.
//
//    Output, bool CH_IS_ALPHANUMERIC is TRUE if the character is alphanumeric.
//
{
  bool value;

  if ( ( 'a' <= ch && ch <= 'z' ) ||
       ( 'A' <= ch && ch <= 'Z' ) ||
       ( '0' <= ch && ch <= '9' ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

bool ch_is_control ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_CONTROL is TRUE if a character is a control character.
//
//  Discussion:
//
//    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to be tested.
//
//    Output, bool CH_IS_CONTROL, TRUE if the character is a control
//    character, and FALSE otherwise.
//
{
  bool value;

  if ( ch <= 31 || 127 <= ch )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
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
//****************************************************************************80

bool ch_is_format_code ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_FORMAT_CODE returns TRUE if a character is a FORTRAN format code.
//
//  Discussion:
//
//    The format codes accepted here are not the only legal format
//    codes in FORTRAN90.  However, they are more than sufficient
//    for my needs!
//
//  Table:
//
//    A  Character
//    B  Binary digits
//    D  Real number, exponential representation
//    E  Real number, exponential representation
//    F  Real number, fixed point
//    G  General format
//    I  Integer
//    L  Logical variable
//    O  Octal digits
//    Z  Hexadecimal digits
//    *  Free format
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_FORMAT_CODE, is TRUE if C is a FORTRAN format code.
//
{
       if ( ch_eqi ( c, 'A' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'B' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'D' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'E' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'F' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'G' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'I' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'L' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'O' ) )
  {
    return true;
  }
  else if ( ch_eqi ( c, 'Z' ) )
  {
    return true;
  }
  else if ( c == '*' )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80

bool ch_is_lower ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_LOWER is TRUE if C is a lowercase alphabetic character.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, a character to check.
//
//    Output, bool CH_IS_LOWER is TRUE if C is a lowercase alphabetic character.
//
{
  bool value;

  if ( ( 'a' <= c && c <= 'z' ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

bool ch_is_printable ( unsigned char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_PRINTABLE determines if a character is printable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned char CH, a character to check.
//
//    Output, bool CH_IS_PRINTABLE is TRUE if the character is printable.
//
{
  bool value;

  if ( 32 <= ch && ch <= 127 )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool ch_is_space ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_SPACE is TRUE if a character represents "white space".
//
//  Discussion:
//
//    A white space character is a space, a form feed, a newline, a carriage
//    return, a horizontal tab, or a vertical tab.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_SPACE, is TRUE if C is a whitespace character.
//
{
  bool value;

  if ( c == ' ' )
  {
    value = true;
  }
  else if ( c == '\f' )
  {
    value = true;
  }
  else if ( c == '\n' )
  {
    value = true;
  }
  else if ( c == '\r' )
  {
    value = true;
  }
  else if ( c == '\t' )
  {
    value = true;
  }
  else if ( c == '\v' )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool ch_is_upper ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_UPPER is TRUE if C is an uppercase alphabetic character.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, a character to check.
//
//    Output, bool CH_IS_UPPER is TRUE if C is an uppercase alphabetic
//    character.
//
{
  bool value;

  if ( ( 'A' <= c && c <= 'Z' ) )
  {
    value = true;
  }
  else
  {
    value = false;
  }

  return value;
}
//****************************************************************************80

char ch_low ( char ch )

//****************************************************************************80
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

int ch_pad ( int *char_index, int *null_index, char *s,
  int max_string )

//****************************************************************************80
//
//  Purpose:
//
//    CH_PAD "pads" a character in a string with a blank on either side.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *CHAR_INDEX, the position of the character to be
//    padded.  On output, this is increased by 1.
//
//    Input/output, int *NULL_INDEX, the position of the terminating NULL in
//    the string.  On output, this is increased by 2.
//
//    Input/output, char *S, the string to be manipulated.
//
//    Input, int MAX_STRING, the maximum number of characters that can be stored
//    in S.
//
//    Output, int CH_PAD, is 0 if the operation worked, and 1 otherwise.
//
{
  int i;

  if ( *char_index < 0 ||
      *null_index <= *char_index ||
      max_string - 1 < *char_index )
  {
    return 1;
  }

  if ( max_string - 1 < (*null_index) + 2 )
  {
    return 1;
  }

  for ( i = *null_index + 2; *char_index + 2 < i; i-- )
  {
    s[i] = s[i-2];
  }
  s[*char_index+2] = ' ';
  s[*char_index+1] = s[*char_index];
  s[*char_index] = ' ';

  *char_index = *char_index + 1;
  *null_index = *null_index + 2;

  return 0;
}
//****************************************************************************80

char ch_read ( FILE *filein )

//****************************************************************************80
//
//  Purpose:
//
//    CH_READ reads one character from a binary file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, FILE *FILEIN, a pointer to the file.
//
//    Output, char CH_READ, the character that was read.
//
{
  char c;

  c = ( char ) fgetc ( filein );

  return c;
}
/******************************************************************************/

int ch_roman_to_i4 ( char ch )

/******************************************************************************/
/*
  Purpose:

    CH_ROMAN_TO_I4 returns the integer value of a single Roman digit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, char CH, a Roman digit.

    Output, int CH_ROMAN_TO_I4, the value of the Roman
    numeral.  If the Roman numeral was not recognized, 0 is returned.
*/
{
  int value;

  if ( ch == 'M' || ch == 'm' )
  {
    value = 1000;
  }
  else if ( ch == 'D' || ch == 'd' )
  {
    value = 500;
  }
  else if ( ch == 'C' || ch == 'c' )
  {
    value = 100;
  }
  else if ( ch == 'L' || ch == 'l' )
  {
    value = 50;
  }
  else if ( ch == 'X' || ch == 'x' )
  {
    value = 10;
  }
  else if ( ch == 'V' || ch == 'v' )
  {
    value = 5;
  }
  else if ( ch == 'I' || ch == 'i' || ch == 'J' || ch == 'j' )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
//****************************************************************************80

char ch_scrabble ( int tile )

//****************************************************************************80
//
//  Purpose:
//
//    CH_SCRABBLE returns the character on a given Scrabble tile.
//
//  Discussion:
//
//    The tiles are numbered 1 to 100, and are labeled 'A' through 'Z',
//    plus two blanks.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TILE, the index of the desired Scrabble tile;
//    1 <= TILE <= 100.
//
//    Output, char CH_SCRABBLE, the character on the given tile.
//
{
  char scrabble[100] = {
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'B',
    'B', 'C', 'C', 'D', 'D', 'D', 'D', 'E', 'E', 'E',
    'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'E', 'F',
    'F', 'G', 'G', 'G', 'H', 'H', 'I', 'I', 'I', 'I',
    'I', 'I', 'I', 'I', 'I', 'J', 'K', 'L', 'L', 'L',
    'L', 'M', 'M', 'N', 'N', 'N', 'N', 'N', 'N', 'O',
    'O', 'O', 'O', 'O', 'O', 'O', 'O', 'P', 'P', 'Q',
    'R', 'R', 'R', 'R', 'R', 'R', 'S', 'S', 'S', 'S',
    'T', 'T', 'T', 'T', 'T', 'T', 'U', 'U', 'U', 'U',
    'V', 'V', 'W', 'W', 'X', 'X', 'Y', 'Z', ' ', ' ' };
  int value;

  if ( 1 <= tile && tile <= 100 )
  {
    value = scrabble[tile-1];
  }
  else
  {
    value = '?';
  }

  return value;
}
//****************************************************************************80

int ch_scrabble_frequency ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_SCRABBLE_FREQUENCY returns the Scrabble frequency of a character.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character.
//
//    Output, int CH_SCRABBLE_FREQUENCY, the frequency of
//    the character.
//
{
  int frequency[27] = {
     9,  2,  2,  4, 12, 
     2,  3,  2,  9,  1, 
     1,  4,  2,  6,  8, 
     2,  1,  6,  4,  6, 
     4,  2,  2,  1,  2, 
     1,  2 };
  int ic;
  int value;
//
//  Convert character to a Scrabble character index.
//
  ic = ch_to_scrabble ( ch );

  if ( 1 <= ic && ic <= 27 )
  {
    value = frequency[ic-1];
  }
  else
  {
    value = 0;
  }

  return value;
}
//****************************************************************************80

int ch_scrabble_points ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_SCRABBLE_POINTS returns the Scrabble point value of a character.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character.
//
//    Output, int CH_SCRABBLE_POINTS, the point value of
//    the character.
//
{
  int ic;
  int points[27] = {
     1,  3,  3,  2,  1, 
     4,  2,  4,  1,  8, 
     5,  1,  3,  1,  1, 
     3, 10,  1,  1,  1, 
     1,  4,  4,  8,  4, 
    10,  0 };
  int value;
//
//  Convert character to a Scrabble character index.
//
  ic = ch_to_scrabble ( ch );

  if ( 1 <= ic && ic <= 27 )
  {
    value = points[ic-1];
  }
  else
  {
    value = 0;
  }

  return value;
}
//****************************************************************************80

char ch_scrabble_select ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CH_SCRABBLE_SELECT selects a character with the Scrabble probability.
//
//  Discussion:
//
//    There are 100 Scrabble tiles, including two blanks.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
//    Output, char CH_SCRABBLE_SELECT, the character on a randomly
//    chosen Scrabble tile.
//
{
  int tile;
  char value;
//
//  Choose a tile between 1 and 100.
//
  tile = i4_uniform ( 1, 100, seed );
//
//  Retrieve the character on that tile.
//
  value = ch_scrabble ( tile );

  return value;
}
//****************************************************************************80

void ch_swap ( char *ch1, char *ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_SWAP swaps two characters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *CH1, *CH2.  On output, the values have been
//    interchanged.
//
{
  char ch3;

   ch3 = *ch1;
  *ch1 = *ch2;
  *ch2 =  ch3;

  return;
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

int ch_to_digit_bin ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT_BIN returns the integer value of a binary digit.
//
//  Discussion:
//
//    This routine handles other traditional binary pairs of "digits"
//    besides '0' and '1'.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    'T'    1
//    'F'    0
//    'Y'    1
//    'N'    0
//    '+'    1
//    '-'    0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the binary digit.
//
//    Output, int CH_TO_DIGIT_BIN, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( c == '0' ||
      c == 'F' ||
      c == 'f' ||
      c == '-' ||
      c == 'N' ||
      c == 'n' )
  {
    digit = 0;
  }
  else if ( c == '1' ||
            c == 'T' ||
            c == 't' ||
            c == '+' ||
            c == 'Y' ||
            c == 'y' )
  {
    digit = 1;
   }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

int ch_to_digit_oct ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT_OCT returns the integer value of an octal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the octal digit, '0' through '7'.
//
//    Output, int CH_TO_DIGIT_OCT, the corresponding integer value, or
//    -1 if C was illegal.
//
{
  int i;

  if ( '0' <= c && c <= '7' )
  {
    i = ( int ) ( c - '0' );
  }
  else if ( c == ' ' )
  {
    i = 0;
   }
  else
  {
    i = -1;
  }

  return i;
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

int ch_to_scrabble ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_SCRABBLE returns the Scrabble index of a character.
//
//  Discussion:
//
//    'A' through 'Z' have indices 1 through 26, and blank is index 27.
//    Case is ignored.  All other characters return index -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character.
//
//    Output, int CH_TO_SCRABBLE, the Scrabble index of
//    the character.
//
{
  int ic;
  int value;

  if ( ch == ' ' )
  {
    value = 27;
    return value;
  }

  ch = ch_cap ( ch );
  ic = a_to_i4 ( ch );

  if ( 1 <= ic && ic <= 26 )
  {
    value = ic;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

char ch_uniform ( char clo, char chi, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CH_UNIFORM returns a random character in a given range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CLO, CHI, the minimum and maximum acceptable characters.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, char CH_UNIFORM, the randomly chosen character.
//
{
  char c;
  double d;

  d = r8_uniform_01 ( seed );

  c = ( char ) ( ( 1.0 - d ) * ( double ) clo + d * ( double ) chi );

  return c;
}
//****************************************************************************80

int ch_write ( FILE *fileout, char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_WRITE writes one character to a binary file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, FILE *FILEOUT, a pointer to the file.
//
//    Input, char C, the character to be written to the file.
//
{
  fputc ( c, fileout );

  return 1;
}
//****************************************************************************80

void charstar_adjustl ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    CHARSTAR_ADJUSTL flushes a CHAR* string left.
//
//  Discussion:
//
//    Both blanks and tabs are treated as "white space".
//
//    This routine is similar to the FORTRAN90 ADJUSTL routine.
//
//  Example:
//
//    Input             Output
//
//    '     Hello'      'Hello'
//    ' Hi there!  '    'Hi there!'
//    'Fred  '          'Fred'
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
//    Input/output, char *S, the string to be adjusted.
//
{
  int i;
  int length;
  int nonb;
  char TAB = 9;
/*
  Check the length of the string to the last nonblank.
  If nonpositive, return.
*/
  length = charstar_len_trim ( s );

  if ( length <= 0 )
  {
    return;
  }
/*
  Find NONB, the location of the first nonblank, nontab.
*/
  nonb = 0;

  for ( i = 0; i < length; i++ )
  {
    if ( s[i] != ' ' && s[i] != TAB )
    {
      nonb = i;
      break;
    }
  }

  if ( 0 < nonb )
  {
    for ( i = nonb; i < length; i++ )
    {
      s[i-nonb] = s[i];
    }

    s[length-nonb] = '\0';
  }
  return;
}
//****************************************************************************80

char *charstar_cat ( char *s1, char *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    CHARSTAR_CAT concatenates two CHAR*'s to make a third.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S1, the "prefix" string.
//
//    Input, char *S2, the "postfix" string.
//
//    Output, char *CHARSTAR_CAT, the string made by
//    concatenating S1 and S2, ignoring any trailing blanks.
//
{
  int i;
  int l1;
  int l2;
  char *s3;

  l1 = charstar_len_trim ( s1 );
  l2 = charstar_len_trim ( s2 );

  if ( l1 == 0 && l2 == 0 )
  {
    s3 = NULL;
    return s3;
  }

  s3 = new char ( l1 + l2 + 1 );

  for ( i = 0; i < l1; i++ )
  {
    s3[i] = s1[i];
  }

  for ( i = 0; i < l2; i++ )
  {
    s3[l1+i] = s2[i];
  }

  s3[l1+l2] = '\0';

  return s3;
}
//****************************************************************************80

bool charstar_eqi ( char *s1, char *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    CHARSTAR_EQI reports whether two CHAR*'s are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S1, char *S2, pointers to two strings.
//
//    Output, bool CHARSTAR_EQI, is true if the strings are equal.
//
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  if ( nchar1 < nchar2 )
  {
    nchar = nchar1;
  }
  else
  {
    nchar = nchar2;
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
  if ( nchar < nchar1 )
  {
    for ( i = nchar; i < nchar1; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return false;
      }
    }
  }
  else if ( nchar < nchar2 )
  {
    for ( i = nchar; i < nchar2; i++ )
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

int charstar_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    CHARSTAR_LEN_TRIM returns the length of a CHAR* to the last nonblank.
//
//  Discussion:
//
//    This function used to be called S_LEN_TRIM.  However, it seems preferable
//    to use the STRING class rather than CHAR*.
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
//    Input, char *S, a pointer to a string.
//
//    Output, int CHARSTAR_LEN_TRIM, the length of the string to the last nonblank.
//    If CHARSTAR_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

char digit_bin_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_BIN_TO_CH returns the character representation of a binary digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer, between 0 and 1.
//
//    Output, char DIGIT_BIN_TO_CH, the character representation of the integer.
//
{
  char c;

  if ( i == 0 )
  {
    c = '0';
  }
  else if ( i == 1 )
  {
    c = '1';
  }
  else
  {
    c = '*';
  }

  return c;
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

char digit_oct_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_OCT_TO_CH returns the character representation of an octal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer, between 0 and 7.
//
//    Output, char DIGIT_OCT_TO_CH, the character representation of the integer.
//
{
  char c;

  if ( 0 <= i && i <= 7 )
  {
    c = i + '0';
  }
  else
  {
    c = '*';
  }

  return c;
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
//    Output, char DIGIT_TO_CH, the appropriate character '0'
//    through '9' or '*'.
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

unsigned int getbits ( unsigned int x, int p, int n )

//****************************************************************************80
//
//  Purpose:
//
//    GETBITS returns N bits from an unsigned int X, beginning at position P.
//
//  Discussion:
//
//    Bits are numbered from right to left, 0 to 31:
//
//       3         2         1
//      10987654321098765432109876543210
//
//    The bits to be extracted begin at position P, and extend TO THE RIGHT.
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
//  Reference:
//
//    Brian Kernighan and Dennis Ritchie,
//    The C Programming Language,
//    Second Edition,
//    Prentice Hall, 1988, page 49.
//
//  Parameters:
//
//    Input, unsigned int X, (or any item of the same size), from which the
//    bits are to be extracted.
//
//    Input, int P, the position in X at which the extraction is to begin.
//    Legal values of P are between 31 (leftmost bit) and 0 (rightmost).
//
//    Input, int N, the number of bits to extract.  N must be at least 1, and
//    P+1-N must be nonnegative.
//
//    Output, int GETBITS, the extracted bits.
//
{
  return ( x >> (p+1-n) ) & ~(~0 << n );
}
//****************************************************************************80

int hex_digit_to_i4 ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_DIGIT_TO_I4 converts a hexadecimal digit to an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the hexadecimal digit, '0'
//    through '9', or 'A' through 'F', or also 'a' through 'f'
//    are allowed.
//
//    Output, int HEX_DIGIT_TO_I4, the corresponding integer,
//    or -1 if C was illegal.
//
{
  int i;

  if ( '0' <= c && c <= '9' )
  {
    i = ( int ) ( c - '0' );
  }
  else if ( 'A' <= c && c <= 'F' )
  {
    i = 10 + ( int ) ( c - 'A' );
  }
  else if ( 'a' <= c && c <= 'f' )
  {
    i = 10 + ( int ) ( c - 'a' );
  }
  else if ( c == ' ' )
  {
    i = 0;
  }
  else
  {
    i = -1;
  }

  return i;
}
//****************************************************************************80

string hex_to_binary_digits ( char hex_digit )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_TO_BINARY_DIGITS converts a hexadecimal digit to 4 binary digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char HEX_DIGIT, the hexadecimal digit.
//
//    Output, character ( len = 4 ) BINARY_DIGITS, the binary digits.
//
{
  string binary_digits;

  if ( hex_digit == '0' )
  {
    binary_digits = "0000";
  }
  else if ( hex_digit == '1' )
  {
    binary_digits = "0001";
  }
  else if ( hex_digit == '2' )
  {
    binary_digits = "0010";
  }
  else if ( hex_digit == '3' )
  {
    binary_digits = "0011";
  }
  else if ( hex_digit == '4' )
  {
    binary_digits = "0100";
  }
  else if ( hex_digit == '5' )
  {
    binary_digits = "0101";
  }
  else if ( hex_digit == '6' )
  {
    binary_digits = "0110";
  }
  else if ( hex_digit == '7' )
  {
    binary_digits = "0111";
  }
  else if ( hex_digit == '8' )
  {
    binary_digits = "1000";
  }
  else if ( hex_digit == '9' )
  {
    binary_digits = "1001";
  }
  else if ( hex_digit == 'A' || hex_digit == 'a' )
  {
    binary_digits = "1010";
  }
  else if ( hex_digit == 'B' || hex_digit == 'b' )
  {
    binary_digits = "1011";
  }
  else if ( hex_digit == 'C' || hex_digit == 'c' )
  {
    binary_digits = "1100";
  }
  else if ( hex_digit == 'D' || hex_digit == 'd' )
  {
    binary_digits = "1101";
  }
  else if ( hex_digit == 'E' || hex_digit == 'e' )
  {
    binary_digits = "1110";
  }
  else if ( hex_digit == 'F' || hex_digit == 'f' )
  {
    binary_digits = "1111";
  }
  else
  {
    binary_digits = "    ";
  }

  return binary_digits;
}
//****************************************************************************80

int hex_to_i4 ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    HEX_TO_I4 converts a hexadecimal string to its integer value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string of hexadecimal digits.
//
//    Output, int HEX_TO_I4, the corresponding integer value.
//
{
  int first;
  int i4;
  int idig;
  int isgn;
  int j;
  int s_length;

  s_length = s.length ( );
//
//  Determine if there is a plus or minus sign.
//
  isgn = 1;

  first = s_length + 1;

  for ( j = 0; j < s_length; j++ )
  {
    if ( s[j] == '-' )
    {
      isgn = - 1;
    }
    else if ( s[j] == '+' )
    {
      isgn = + 1;
    }
    else if ( s[j] != ' ' )
    {
      first = j;
      break;
    }
  }
//
//  Read the numeric portion of the string.
//
  i4 = 0;

  for ( j = first; j < s_length; j++ )
  {
    idig = hex_digit_to_i4 ( s[j] );
    i4 = i4 * 16 + idig;
  }

  i4 = isgn * i4;

  return i4;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4, usually the largest legal signed int.
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
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80

int i4_input ( string prompt, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    I4_INPUT prints a prompt string and reads an I4 from the user.
//
//  Discussion:
//
//    If the input line starts with a comment character ('#') or is
//    blank, the routine ignores that line, and tries to read the next one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PROMPT, the prompt string.
//
//    Output, bool &ERROR, an error flag, which is true if an error occurred.
//
//    Output, integer I4_INPUT, the value input by the user.
//
{
  int last;
  char line[80];
  int value;

  error = false;
  value = i4_huge ( );
//
//  Write the prompt.
//
  cout << "\n";
  cout << prompt << "\n";

  for ( ; ; )
  {
    cin.getline ( line, sizeof ( line ) );
//
//  If the line begins with a comment character, go back and read the next line.
//
    if ( line[0] == '#' )
    {
      continue;
    }

    if ( charstar_len_trim ( line ) == 0 )
    {
      continue;
    }
//
//  Extract integer information from the string.
//
    value = s_to_i4 ( line, last, error );

    if ( error )
    {
      value = i4_huge ( );
      return value;
    }
    break;
  }

  return value;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the whole part of the logarithm base 10 of an I4.
//
//  Discussion:
//
//    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
//    (except for I = 0).
//
//    The number of decimal digits in I is I4_LOG_10(I) + 1.
//
//  Example:
//
//        I    I4_LOG_10(I)
//
//        0     0
//        1     0
//        2     0
//
//        9     0
//       10     1
//       11     1
//
//       99     1
//      100     2
//      101     2
//
//      999     2
//     1000     3
//     1001     3
//
//     9999     3
//    10000     4
//    10001     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer.
//
//    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
//
{
  int ten_pow;
  int value;

  i = abs ( i );

  ten_pow = 10;
  value = 0;

  while ( ten_pow <= i )
  {
    ten_pow = ten_pow * 10;
    value = value + 1;
  }

  return value;
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
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
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
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
//****************************************************************************80

char i4_to_a ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_A returns the I-th alphabetic character.
//
//  Example:
//
//    I  I4_TO_A
//
//   -8  ' '
//    0  ' '
//    1  'A'
//    2  'B'
//   ..
//   26  'Z'
//   27  'a'
//   52  'z'
//   53  ' '
//   99  ' '
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the index of the letter to be returned.
//    0 is a space;
//    1 through 26 requests 'A' through 'Z', (ASCII 65:90);
//    27 through 52 requests 'a' through 'z', (ASCII 97:122);
//
//    Output, char I4_TO_A, the requested alphabetic letter.
//
{
  char value;

  if ( i <= 0 )
  {
    value = ' ';
  }
  else if ( 1 <= i && i <= 26 )
  {
    value = 'A' + i - 1;
  }
  else if ( 27 <= i && i <= 52 )
  {
    value = 'a' + i - 27;
  }
  else
  {
    value = ' ';
  }
  return value;
}
//****************************************************************************80

char i4_to_amino_code ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_AMINO_CODE converts an integer to an amino code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl Branden, John Tooze,
//    Introduction to Protein Structure,
//    Garland Publishing, 1991.
//
//  Parameters:
//
//    Input, int I, the index of an amino acid, between 1 and 23.
//
//    Output, char I4_TO_AMINO_CODE, the one letter code for an amino acid.
//
{
# define N 23

  char c;
  static char ch_table[N] = {
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
    'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
    'X', 'Y', 'Z' };

  if ( 1 <= i && i <= N )
  {
    c = ch_table[i-1];
  }
  else
  {
    c = '?';
  }

  return c;
# undef N
}
//****************************************************************************80

char i4_to_hex_digit ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HEX_DIGIT converts a (small) I4 to a hexadecimal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer, between 0 and 15.
//
//    Output, char DI4_TO_HEX_DIGIT, the hexadecimal digit corresponding
//    to the integer.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else if ( 10 <= i && i <= 15 )
  {
    c = 'a' + ( i - 10 );
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

char i4_to_isbn ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_ISBN converts an I4 to an ISBN digit.
//
//  Discussion:
//
//    Only the integers 0 through 10 can be input.  The representation
//    of 10 is 'X'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input, int I, an integer between 0 and 10.
//
//    Output, char I4_TO_ISBN, the ISBN character code of the integer.
//    If I is illegal, then I4_TO_ISBN is set to '?'.
//
{
       if ( i == 0 )
  {
    return '0';
  }
  else if ( i == 1 )
  {
    return '1';
  }
  else if ( i == 2 )
  {
    return '2';
  }
  else if ( i == 3 )
  {
    return '3';
  }
  else if ( i == 4 )
  {
    return '4';
  }
  else if ( i == 5 )
  {
    return '5';
  }
  else if ( i == 6 )
  {
    return '6';
  }
  else if ( i == 7 )
  {
    return '7';
  }
  else if ( i == 8 )
  {
    return '8';
  }
  else if ( i == 9 )
  {
    return '9';
  }
  else if ( i == 10 )
  {
    return 'X';
  }
  else
  {
    return '?';
  }
}
//****************************************************************************80

string i4_to_month_abb ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_MONTH_ABB returns an abbreviated month name.
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
//    Input, int I, the number of the desired month.
//
//    Output, string I4_TO_MONTH_ABB, a 3 character abbreviation for
//    the month, such as 'Jan', 'Feb', 'Mar', and so on.
//
{
  string s;
  static string month_list[13] =
  {
    "???",
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  };

  if ( 1 <= i && i <= 12 )
  {
    s = month_list[i];
  }
  else
  {
    s = month_list[0];
  }
  return s;
}
//****************************************************************************80

string i4_to_month_name ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_MONTH_NAME returns a month name.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number of the desired month.
//    1 <= I <= 12.
//
//    Output, string I4_TO_MONTH_NAME, the name of the month.
//
{
  string value;

  if ( i == 1 )
  {
    value = "January";
  }
  else if ( i == 2 )
  {
    value = "February";
  }
  else if ( i == 3 )
  {
    value = "March";
  }
  else if ( i == 4 )
  {
    value = "April";
  }
  else if ( i == 5 )
  {
    value = "May";
  }
  else if ( i == 6 )
  {
    value = "June";
  }
  else if ( i == 7 )
  {
    value = "July";
  }
  else if ( i == 8 )
  {
    value = "August";
  }
  else if ( i == 9 )
  {
    value = "September";
  }
  else if ( i == 10 )
  {
    value = "October";
  }
  else if ( i == 11 )
  {
    value = "November";
  }
  else if ( i == 12 )
  {
    value = "December";
  }
  else
  {
    value = "???";
  }
  return value;
}
//****************************************************************************80

string i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, string I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  string s;
  char s_char[80];
  static double ten = 10.0;

  length = i4_log_10 ( i );

  ten_power = ( int ) ( pow ( ten, length ) );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  if ( i == 0 )
  {
    s_char[0] = '0';
    s_char[1] = '\0';
    s = string ( s_char );
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s_char[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s_char[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s_char[j] = '\0';
  j = j + 1;

  s = string ( s_char );

  return s;
}
//****************************************************************************80

string i4_to_s0 ( int i, int digits )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S0 converts an I4 to a string with leading zeros.
//
//  Example:
//
//         I  DIGITS       S
//
//         1       3     001
//        -1       3     -01
//         0       3     000
//      1952       4    1952
//    123456       6  123456
//   1234567       6  ******
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Input, int DIGITS, the number of positions available.
//
//    Output, string I4_TO_S0, the representation of the integer.
//
{
  int j;
  int length;
  int length_max = 80;
  string s;
  char s_char[81];
//
//  Determine LEGNTH, the natural "length" of the number.
//
  length = i4_log_10 ( i ) + 1;

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  If DIGITS is too large, reset it.
//
  digits = i4_min ( digits, length_max );
//
//  If I is too large, return stars.
//
  if ( digits < length )
  {
    for ( j = 0; j < digits; j++ )
    {
      s_char[j] = '*';
    }
    s_char[digits] = '\0';
    s = string ( s_char );
    return s;
  }
//
//  Set the digits to zero.
//
  for ( j = 0; j < digits; j++ )
  {
    s_char[j] = '0';
  }
  s_char[digits] = '\0';
//
//  If I = 0, we can return now.
//
  if ( i == 0 )
  {
    s = string ( s_char );
    return s;
  }
//
//  If negative, take care of the sign.
//
  if ( i < 0 )
  {
    s_char[0] = '-';
    i = abs ( i );
  }
//
//  Mod away.
//
  j = digits - 1;
  while ( i != 0 )
  {
    s_char[j] = digit_to_ch ( i % 10 );
    i = i / 10;
    j = j - 1;
  }
  s = string ( s_char );

  return s;
}
//****************************************************************************80

string i4_to_string ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << i4;

  value = fred.str ( );

  return value;
}
//****************************************************************************80

string i4_to_unary ( int i4 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_UNARY produces the "base 1" representation of an I4.
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
//    Input, int I4, an integer to be represented.
//
//    Output, string S, the unary representation.
//
{
  int i;
  string s;
  char *s_char;

  if ( i4 < 0 )
  {
    s_char = new char[-i4+2];
    s_char[0] = '-';
    for ( i = 1; i <= -i4; i++ )
    {
      s_char[i] = '1';
    }
    s_char[-i4+1] = '\0';
  }
  else if ( i4 == 0 )
  {
    s_char = new char[2];
    s_char[0] = '0';
    s_char[1] = '\0';
  }
  else if ( 0 < i4 )
  {
    s_char = new char[i4+1];
    for ( i = 0; i < i4; i++ )
    {
      s_char[i] = '1';
    }
    s_char[i4] = '\0';
  }

  s = string ( s_char );

  delete [] s_char;

  return s;
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

int *i4vec_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
//
//  Discussion:
//
//    An I4VEC is a vector of integer values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR(N), the initialized array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }

  return a;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of integer values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i + 1 << "  "
         << setw(8) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void int_to_byte ( unsigned int ival, unsigned char *bvec )

//****************************************************************************80
//
//  Purpose:
//
//    INT_TO_BYTE converts an unsigned integer into 4 bytes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned int IVAL, is the integer to be converted.
//
//    Output, unsigned char *BVEC, is a pointer to a character string.
//    The contents of BVEC through BVEC+3 will be the bytes of IVAL,
//    from low order to high.
//
{
  int i;

  for ( i = 0; i < 4; i++ )
  {
    *bvec = ( ival >> (3-i)*8 );
    bvec = bvec + 1;
  }
  return;
}
//****************************************************************************80

int isbn_to_i4 ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    ISBN_TO_I4 converts an ISBN character into an I4.
//
//  Discussion:
//
//    The characters '0' through '9' stand for themselves, but
//    the character 'X' or 'x' stands for 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Book Industry Study Group,
//    The Evolution in Product Identification:
//    Sunrise 2005 and the ISBN-13,
//    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
//
//  Parameters:
//
//    Input, char C, the ISBN character code to be converted.
//
//    Output, int ISBN_TO_I4, the numeric value of the character
//    code, between 0 and 10.  This value is returned as -1 if C is
//    not a valid character code.
//
{
  int value;

  if ( '0' <= c && c <= '9' )
  {
    value = c - '0';
  }
  else if ( c == 'X' || c == 'x' )
  {
    value = 10;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

bool perm_check ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = 1; seek <= n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//****************************************************************************80

int *perm_uniform ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM[N], a permutation of (1,, 1, ..., N).
//
{
  int i;
  int j;
  int *p;

  p = new int[n];

  for ( i = 1; i <= n; i++ )
  {
    p[i-1] = i;
  }

  for ( i = 1; i <= n; i++ )
  {
    j = i4_uniform ( i, n, seed );
    i4_swap ( &p[i-1], &p[j-1] );
  }

  return p;
}
//****************************************************************************80

void print_sizes ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_SIZES reports the size in bytes of various data types.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
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
  cout << "\n";
  cout << "PRINT_SIZES: Report data type sizes.\n";
  cout << "\n";
  cout << "  (Min/Max data currently unavailable...\n";
  cout << "\n";
  cout << "  Type                   Size     Min     Max\n";
  cout << "\n";
  cout << "  bool                   " << sizeof ( bool ) << "\n";
  cout << "  char                   " << sizeof ( char ) << "\n";
  cout << "  unsigned char          " << sizeof ( char )
    << " \n";
  cout << "  short int              " << sizeof ( short int ) << "\n";
  cout << "  unsigned short int     " << sizeof ( unsigned short int )
       << " \n";
  cout << "  int                    " << sizeof ( int ) << "\n";
  cout << "  unsigned int           " << sizeof ( unsigned int )
       << " \n";
  cout << "  long int               " << sizeof ( long int ) << "\n";
  cout << "  unsigned long int      " << sizeof ( unsigned long int )
       << " \n";

  cout << "  long long int          " << sizeof ( long long int ) << "\n";
  cout << "  unsigned long long int " << sizeof ( unsigned long long int )
       << " \n";

  cout << "  float                  " << sizeof ( float ) << "\n";
  cout << "  double                 " << sizeof ( double ) << "\n";

  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
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

string r4_to_string ( float r4 )

//****************************************************************************80
//
//  Purpose:
//
//    R4_TO_STRING converts an R4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2-13
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float R4, a float.
//
//    Output, string R4_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << r4;

  value = fred.str ( );

  return value;
}
//****************************************************************************80

string r8_to_string ( double r8 )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_STRING converts an R8 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 january 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R8, a double.
//
//    Output, string R8_TO_STRING, the string.
//
{
  ostringstream fred;
  string value;

  fred << r8;

  value = fred.str ( );

  return value;
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
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

float reverse_bytes_float ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    REVERSE_BYTES_FLOAT reverses the four bytes in a float.
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
//    Input, float X, a float whose bytes are to be reversed.
//
//    Output, float REVERSE_BYTES_FLOAT, a float with bytes in reverse order
//    from those in x.
//
{
  char c;
  union
  {
    float yfloat;
    char ychar[4];
  } y;

  y.yfloat = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[3];
  y.ychar[3] = c;

  c = y.ychar[1];
  y.ychar[1] = y.ychar[2];
  y.ychar[2] = c;

  return ( y.yfloat );
}
//****************************************************************************80

int reverse_bytes_int ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    REVERSE_BYTES_INT reverses the four bytes in an int.
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
//    Input, int X, an int whose bytes are to be reversed.
//
//    Output, int REVERSE_BYTES_INT, an int with bytes in reverse order
//    from those in X.
//
{
  char c;
  union
  {
    int yint;
    char ychar[4];
  } y;

  y.yint = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[3];
  y.ychar[3] = c;

  c = y.ychar[1];
  y.ychar[1] = y.ychar[2];
  y.ychar[2] = c;

  return ( y.yint );
}
//****************************************************************************80

string s_adjustl ( string s1 )

//****************************************************************************80
//
//  Purpose:
//
//    S_ADJUSTL flushes a string left.
//
//  Discussion:
//
//    Both blanks and tabs are treated as "white space".
//
//    This routine is similar to the FORTRAN90 ADJUSTL routine.
//
//  Example:
//
//    Input             Output
//
//    '     Hello'      'Hello'
//    ' Hi there!  '    'Hi there!'
//    'Fred  '          'Fred'
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
//    Input, string S1, the string to be adjusted.
//
//    Output, string S_ADJUSTL, the adjusted string.
//
{
  int i;
  int s2_length;
  string s2;
  int nonb;
  char TAB = 9;

  s2 = s1;
//
//  Check the length of the string to the last nonblank.
//  If nonpositive, return.
//
  s2_length = s2.length ( );

  if ( s2_length <= 0 )
  {
    return s2;
  }
//
//  Find NONB, the location of the first nonblank, nontab.
//
  nonb = 0;

  for ( i = 0; i < s2_length; i++ )
  {
    if ( s1[i] != ' ' && s1[i] != TAB )
    {
      nonb = i;
      break;
    }
  }

  if ( 0 < nonb )
  {
    for ( i = 0; i < s2_length - nonb; i++ )
    {
      s2[i] = s1[i+nonb];
    }
    for ( i = s2_length - nonb; i < s2_length; i++ )
    {
      s2[i] = ' ';
    }

  }
  return s2;
}
//****************************************************************************80

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
{
  int i;
  int n1;
  int n2;

  n1 = s1.length ( );
  n2 = s2.length ( );

  if ( n1 < n2 )
  {
    return false;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

void s_behead_substring ( char *s, char *sub )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEHEAD_SUBSTRING "beheads" a string, removing a given substring.
//
//  Discussion:
//
//    Initial blanks in the string are removed first.
//
//    Then, if the initial part of the string matches the substring,
//    that part is removed and the remainder shifted left.
//
//    Initial blanks in the substring are NOT ignored.
//
//    Capitalization is ignored.
//
//    If the substring is equal to the string, then the resultant
//    string is returned as a single blank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *S, the string to be transformed.
//
//    Input, char *SUB, the substring to be removed from
//    the beginning of the string.
//
{
  int i;
  int s_len;
  int sub_len;
//
//  Remove leading blanks from the string.
//
  charstar_adjustl ( s );
//
//  Get lengths.
//
  s_len = charstar_len_trim ( s );
  sub_len = charstar_len_trim ( sub );

  if ( s_len < sub_len )
  {
    return;
  }
//
//  Does the string match the substring?
//
  for ( i = 0; i < sub_len; i++ )
  {
    if ( !ch_eqi ( *(s+i), *(sub+i) ) )
    {
      return;
    }
  }
//
//  Blank out the substring.
//
  if ( sub_len < s_len )
  {
    for ( i = 0; i < sub_len; i++ )
    {
      *(s+i) = ' ';
    }
    charstar_adjustl ( s );
  }
  else
  {
    strcpy ( s, " " );
  }

  return;
}
//****************************************************************************80

void s_blank_delete ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_BLANK_DELETE removes blanks and left justifies the remainder.
//
//  Discussion:
//
//    All TAB characters are also removed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *S, the string to be transformed.
//
{
  char *get;
  char *put;
  char TAB = 9;

  put = s;
  get = s;

  while ( *get != '\0' )
  {
    if ( *get != ' ' && *get != TAB )
    {
      *put = *get;
      put = put + 1;
    }
    get = get + 1;
  }

  *put = *get;

  return;
}
//****************************************************************************80

string s_blanks_delete ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_BLANKS_DELETE replaces consecutive blanks by one blank.
//
//  Discussion:
//
//    The remaining characters are left justified and right padded with blanks.
//    TAB characters are converted to spaces.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be transformed.
//
//    Output, string S_BLANKS_DELETE, the transformed string.
//
{
  bool blank;
  int get;
  int put;
  int s_length;
  char *s2;
  string s3;

  s_length = s.length ( );
  s2 = new char[s_length+1];
  s2[s_length] = '\0';

  blank = true;
  put = 0;

  for ( get = 0; get < s_length; get++ )
  {
    if ( s[get] != ' ' )
    {
      s2[put] = s[get];
      put = put + 1;
      blank = false;
    }
    else if ( !blank )
    {
      s2[put] = s[get];
      put = put + 1;
      blank = true;
    }
    else
    {
    }
  }
//
//  Suppress a final blank that is not the only character.
//
  if ( 1 < put )
  {
    if ( s2[put-1] == ' ' )
    {
      put = put - 1;
    }
  }
  s2[put] = '\0';
  s3 = string ( s2 );

  return s3;
}
//****************************************************************************80

string s_cap ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_CAP capitalizes all the characters in a string.
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
//    Input, string S, the string to be capitalized.
//
//    Output, string S_CAP, the capitalized string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    s2[i] = ch_cap ( s2[i] );
  }

  return s2;
}
//****************************************************************************80

int s_ch_count ( string s, char ch )

//****************************************************************************80
//
//  Purpose:
//
//    S_CH_COUNT counts occurrences of a particular character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string.
//
//    Input, char CH, the character to be counted.
//
//    Output, int S_CH_COUNT, the number of occurrences.
//
{
  int ch_count;
  int i;
  int s_length;

  s_length = s.length ( );

  ch_count = 0;

  for ( i = 0; i < s_length; i++ )
  {
    if ( s[i] == ch )
    {
      ch_count = ch_count + 1;
    }
  }
  return ch_count;
}
//****************************************************************************80

string s_ch_delete ( string s, char ch )

//****************************************************************************80
//
//  Purpose:
//
//    S_CH_DELETE removes all occurrences of a character from a string.
//
//  Discussion:
//
//    Each time the given character is found in the string, the characters
//    to the right of the string are shifted over one position.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be transformed.
//
//    Input, character CH, the character to be removed.
//
//    Output, string S_CH_DELETE, a copy of the string with the character removed.
//
{
  int ch_num;
  int get;
  int put;
  int s_length;
  string s2;

  s_length = s.length ( );

  ch_num = s_ch_count ( s, ch );

  s2 = string ( s_length - ch_num + 1, ' ' );

  put = 0;

  for ( get = 0; get < s_length; get++ )
  {
    if ( s2[get] != ch )
    {
      s2[put] = s[get];
      put = put + 1;
    }
  }
  return s2;
}
//****************************************************************************80

string s_control_blank ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_CONTROL_BLANK replaces control characters with blanks.
//
//  Discussion:
//
//    A "control character" has ASCII code <= 31 or 127 <= ASCII code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be transformed.
//
//    Output, string S_CONTROL_BLANK, the transformed string.
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    if ( ch_is_control ( s2[i] ) )
    {
      s2[i] = ' ';
    }
  }
  return s2;
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

string s_escape_tex ( string s1 )

//****************************************************************************80
//
//  Purpose:
//
//    S_ESCAPE_TEX de-escapes TeX escape sequences.
//
//  Discussion:
//
//    In particular, every occurrence of the characters '\', '_',
//    '^', '{' and '}' will be replaced by '\\', '\_', '\^',
//    '\{' and '\}'.  A TeX interpreter, on seeing these character
//    strings, is then likely to return the original characters.
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
//    Input, string S1, the string to be de-escaped.
//
//    Output, string S_ESCAPE_TEX, a copy of the string,
//    modified to avoid TeX escapes.
//
{
  char ch;
  int s1_length;
  int s1_pos;
  char *s2;
  int s2_pos;
  string s3;
  int slash_count;

  s1_length = s1.length ( );

  slash_count = 0;
  for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
  {
    ch = s1[s1_pos];

    if ( ch == '\\' ||
         ch == '_' ||
         ch == '^' ||
         ch == '{' ||
         ch == '}' )
    {
      slash_count = slash_count + 1;
    }
  }
  s2 = new char[s1_length + slash_count + 1];
//
//  Now copy S1 into S2.
//
  s1_pos = 0;
  s2_pos = 0;

  for ( s1_pos = 0; s1_pos < s1_length; s1_pos++ )
  {
    ch = s1[s1_pos];

    if ( ch == '\\' ||
         ch == '_' ||
         ch == '^' ||
         ch == '{' ||
         ch == '}' )
    {
      s2[s2_pos] = '\\';
      s2_pos = s2_pos + 1;
    }

    s2[s2_pos] = ch;
    s2_pos = s2_pos + 1;
  }

  s2[s2_pos] = '\0';
  s2_pos = s2_pos + 1;

  s3 = string ( s2 );

  return s3;
}
//****************************************************************************80

int s_first_ch ( string s, char ch )

//****************************************************************************80
//
//  Purpose:
//
//    S_FIRST_CH points to the first occurrence of a character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Input, char CH, a character.
//
//    Output, int S_FIRST_CH, the first occurrence of the
//    character in the string, or -1 if the character does not occur.
//
{
  int i;
  int s_length;
  int value;

  s_length = s.length ( );
  value = - 1;

  for ( i = 0; i < s_length; i++ )
  {
    if ( s[i] == ch )
    {
      value = i;
      return value;
    }
  }

  return value;
}
//****************************************************************************80

char *s_first_nonblank ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_FIRST_NONBLANK points to the first nonblank character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, char *S_FIRST_NONBLANK, a pointer to the first nonblank character
//    in the string, or NULL if the entire string is blank.
//
{
  char *t = NULL;

  while ( *s != '\0' )
  {
    if ( *s != ' ' )
    {
      t = s;
      break;
    }
    s++;
  }

  return t;
}
//****************************************************************************80

void s_inc_c ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_INC_C "increments" the characters in a string.
//
//  Discussion:
//
//    The routine tries to produce the next string, in dictionary order,
//    following the input value of a string.  Digits, spaces, and other
//    nonalphabetic characters are ignored.  Case is respected; in other
//    words, the case of every alphabetic character on input will be the
//    same on output.
//
//    The following error conditions can occur:
//
//      There are no alphabetic characters in the string.  No
//      incrementing is possible.
//
//      All alphabetic characters are equal to 'Z' or 'z'.  In this
//      case, an error value is returned, but the string is also "wrapped
//      around" so that all alphabetic characters are "A" or "a".
//
//    If the word "Tax" were input, the successive outputs would be
//    "Tay", "Taz", "Tba", "Tbb", ...  If the input word "January 4, 1989"
//    were input, the output would be "Januarz 4, 1989".
//
//    This routine could be useful when trying to create a unique file
//    name or variable name at run time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *S, the string whose
//    alphabetic successor is desired.
//
{
  int n;
  char *t;

  n = strlen ( s ) - 1;
  t = s + strlen ( s ) - 1;

  while ( 0 <= n )
  {
    if ( !ch_is_alpha ( *t ) )
    {
    }
    else if ( *t == 'Z' )
    {
      *t = 'A';
    }
    else if ( *t == 'z' )
    {
      *t = 'a';
    }
    else
    {
      *t = *t + 1;
      break;
    }
    t--;
    n--;
  }
  return;
}
//****************************************************************************80

void s_inc_n ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_INC_N increments the digits in a string.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the input string contains no digits, a blank string is returned.
//
//    If a blank string is input, then an error condition results.
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
//    11 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *S, the character string to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = charstar_len_trim ( s );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "S_INC_N - Fatal error!\n";
    cerr << "  Input file name is blank.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens-1; 0 <= i; i-- )
  {
    c = *(s+i);

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      if ( c == '9' )
      {
        c = '0';
        *(s+i) = c;
      }
      else
      {
        c = c + 1;
        *(s+i) = c;
        return;
      }
    }
  }

  if ( change == 0 )
  {
    strcpy ( s, " " );
  }

  return;
}
//****************************************************************************80

string s_last_ch ( string s, char ch )

//****************************************************************************80
//
//  Purpose:
//
//    S_LAST_CH points to the last occurrence of a character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Input, char CH, a character.
//
//    Output, string S_LAST_CH, the substring beginning with the last occurrence
//    of the given character.
//
{
  int position;
  int s_length;
  string t;
  int t_length;

  s_length = s.length ( );
//
//  Find the last occurrence.
//
  for ( position = s_length - 1; 0 <= position; position-- )
  {
    if ( s[position] == ch )
    {
      t_length = s_length - position;
      t = s.substr ( position, t_length );
      return t;
    }
  }

  t.clear ( );

  return t;
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
//    10 October 2014
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
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

string s_low ( string s )

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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be lowercased.
//
//    Output, string S_LOW, the lowercased string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    s2[i] = tolower ( s2[i] );
  }

  return s2;
}
//****************************************************************************80

string s_newline_to_null ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_NEWLINE_TO_NULL replaces carriage returns or newlines by nulls.
//
//  Discussion:
//
//    The function FGETS will read a string containing a line of text read from
//    input.  However, the string will include the linefeed character '/n', or,
//    for a PC-formatted file, the carriage return and linefeed pair '/r' + '/n'.
//
//    It may be desirable that the string not contain these characters.  The
//    easiest way to deal with this is simply to replace the first instance of
//    '/r' or '/n' by a null character, which terminates the string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be modified.
//
//    Output, string S_NEWLINE_TO_NULL, the modified string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
//
//  Handle carriage return;
//
    if ( s[i] == '\r' )
    {
      s2[i] = '\0';
      break;
    }
//
//  Handle linefeed.
//
    if ( s[i] == '\n' )
    {
      s2[i] = '\0';
      break;
    }
  }
  return s2;
}
//****************************************************************************80

string s_nonalpha_delete ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_NONALPHA_DELETE removes nonalphabetic characters from a string.
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
//    Input, string S, the input string.
//
//    Output, string S_NONALPHA_DELETE, the output string.
{
  int i;
  int j;
  int s_length;
  string s2;

  s_length = s.length ( );

  for ( i = 0; i < s_length; i++ )
  {
    if ( ch_is_alpha ( s[i] ) )
    {
      s2[j] = s[i];
      j = j + 1;
    }
  }

  return s2;
}
//****************************************************************************80

string s_replace_ch ( string s, char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_REPLACE_CH replaces all occurrences of one character by another.
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
//    Input, string S, the string.
//
//    Input, char C1, C2, the character to be replaced, and the
//    replacement character.
//
//    Output, string S_REPLACE_CH, the modified string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    if ( s2[i] == c1 )
    {
      s2[i] = c2;
    }
  }
  return s2;
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

bool s_s_subanagram ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_S_SUBANAGRAM determines if S2 is a "subanagram" of S1.
//
//  Discussion:
//
//    S2 is an anagram of S1 if S2 can be formed by permuting the letters
//    of S1
//
//    S2 is an subanagram of S1 if S2 can be formed by selecting SOME of
//    the letters of S1 and permuting them.
//
//    Blanks (trailing or otherwise), punctuation, and capitalization
//    are all significant, so be sure to input exactly the information
//    you want to check.
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
//    Input, string S1, the master string.
//
//    Input, string S2, the second string.
//
//    Output, bool S_S_SUBANAGRAM is TRUE if S2 is a subanagram of S1.
//
{
  int i1;
  int i2;
  int s1_length;
  string s1_sorted;
  int s2_length;
  string s2_sorted;
  bool value;

  value = false;
//
//  Sort both.
//
  s1_sorted = s_sort_a ( s1 );
  s2_sorted = s_sort_a ( s2 );

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  i1 = - 1;

  for ( i2 = 0; i2 < s2_length; i2++ )
  {
    for ( ; ; )
    {
      i1 = i1 + 1;
//
//  Ran out of S1 before finishing.  No match is possible.
//
      if ( s1_length <= i1 )
      {
        return value;
      }
//
//  The current character in S1 is already greater than the character in S2.
//  No match is possible.
//
      if ( s2_sorted[i2] < s1_sorted[i1] )
      {
        return value;
      }
//
//  Found an exact match for current character.  Keep going.
//
      if ( s1_sorted[i1] == s2_sorted[i2] )
      {
        break;
      }
//
//  Didn't find a match, but one might be possible if we increase I1.
//
    }
  }
//
//  We matched every character of S2 with something in S1.
//
  value = true;

  return value;
}
//****************************************************************************80

bool s_s_subanagram_sorted ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_S_SUBANAGRAM_SORTED determines if S2 is a "subanagram" of S1.
//
//  Discussion:
//
//    This function assumes both input strings have already been sorted.
//
//    S2 is an anagram of S1 if S2 can be formed by permuting the letters
//    of S1
//
//    S2 is an subanagram of S1 if S2 can be formed by selecting SOME of
//    the letters of S1 and permuting them.
//
//    Blanks (trailing or otherwise), punctuation, and capitalization
//    are all significant, so be sure to input exactly the information
//    you want to check.
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
//    Input, string S1, the master sorted string.
//
//    Input, string S2, the second sorted string.
//
//    Output, bool S_S_SUBANAGRAM_SORTED is TRUE if S2 is a subanagram of S1.
//
{
  int i1;
  int i2;
  int s1_length;
  int s2_length;
  bool value;

  value = false;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  i1 = - 1;

  for ( i2 = 0; i2 < s2_length; i2++ )
  {
    for ( ; ; )
    {
      i1 = i1 + 1;
//
//  Ran out of S1 before finishing.  No match is possible.
//
      if ( s1_length <= i1 )
      {
        return value;
      }
//
//  The current character in S1 is already greater than the character in S2.
//  No match is possible.
//
      if ( s2[i2] < s1[i1] )
      {
        return value;
      }
//
//  Found an exact match for current character.  Keep going.
//
      if ( s1[i1] == s2[i2] )
      {
        break;
      }
//
//  Didn't find a match, but one might be possible if we increase I1.
//
    }
  }
//
//  We matched every character of S2 with something in S1.
//
  value = true;

  return value;
}
//****************************************************************************80

int s_scrabble_points ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_SCRABBLE_POINTS returns the Scrabble point value of a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string.
//
//    Output, int S_SCRABBLE_POINTS, the point value of
//    the string.
//
{
  int i;
  int s_length;
  int value;

  s_length = s.length ( );

  value = 0;
  for ( i = 0; i < s_length; i++ )
  {
    value = value + ch_scrabble_points ( s[i] );
  }

  return value;
}
//****************************************************************************80

string s_sort_a ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_SORT_A sorts a string into ascending order.
//
//  Discussion:
//
//    The string is assumed to be short, and so a simple bubble sort is used.
//
//    ALL the characters are sorted, including blanks and punctuation.
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
//    Input, string S, the string to be sorted.
//
//    Input, string S2, the sorted string.
//
{
  char c;
  int i;
  int j;
  int k;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length - 1; i++ )
  {
    c = s2[i];
    j = i;

    for ( k = i + 1; k < s_length; k++ )
    {
      if ( s2[k] < s2[j] )
      {
        j = k;
      }
    }

    if ( i != j )
    {
      s2[i] = s2[j];
      s2[j] = c;
    }
  }

  return s2;
}
//****************************************************************************80

char *s_substring ( char *s, int a, int b )

//****************************************************************************80
//
//  Purpose:
//
//    S_SUBSTRING returns a substring of a given string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Input, int A, B, the indices of the first and last character of S to copy.
//    These are 1-based indices!  B should be
//
//    Output, char *S_SUBSTRING, a pointer to the substring.
//
{
  int i;
  int j;
  char *t;

  t = new char[b+2-a];

  j = 0;
  for ( i = a; i <= b; i++ )
  {
    t[j] = s[i-1];
    j = j + 1;
  }
  t[j] = '\0';

  return t;
}
//****************************************************************************80

string s_tab_blank ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_TAB_BLANK replaces each TAB character by a space.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be transformed.
//
//    Output, string S_TAB_BLANK, the transformed string.
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ )
  {
    if ( s2[i] == '\t' )
    {
      s2[i] = ' ';
    }
  }
  return s2;
}
//****************************************************************************80

void s_to_format ( char *s, int *r, char *code, int *w, int *m )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_FORMAT reads a FORTRAN format from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the format.  This routine is limited in its ability to
//    recognize FORTRAN formats.  In particular, we are only expecting
//    a single format specification, and cannot handle extra features
//    such as 'ES' and 'EN' codes, '5X' spacing, and so on.
//
//    Legal input is:
//
//       0 nothing
//       1 blanks
//       2 optional '('
//       3 blanks
//       4 optional repeat factor R
//       5 blanks
//       6 CODE ( 'A', 'B', 'E', 'F', 'G', 'I', 'L', 'O', 'Z', '*' )
//       7 blanks
//       8 width W
//       9 optional decimal point
//      10 optional mantissa M
//      11 blanks
//      12 optional ')'
//      13 blanks
//
//  Example:
//
//    S                 R   CODE   W    M
//
//    'I12              1   I      12   0
//    'E8.0'            1   E       8   0
//    'F10.5'           1   F      10   5
//    '2G14.6'          2   G      14   6
//    '*'               1   *      -1  -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read.
//
//    Output, int *R, the repetition factor, which defaults to 1.
//
//    Output, char *CODE, the format code.
//
//    Output, int *W, the field width.
//
//    Output, int *M, the mantissa width.
//
{
  char c;
  int d;
  bool debug = true;
  int LEFT = 1;
  int paren_sum;
  int pos;
  int RIGHT = -1;
  int s_length;
  int state;

  state = 0;
  paren_sum = 0;
  pos = 0;
  s_length = charstar_len_trim ( s );

  *r = 0;
  *w = 0;
  *code = '?';
  *m = 0;

  while ( pos < s_length )
  {
    c = s[pos];
    pos = pos + 1;
//
//  BLANK character:
//
    if ( c == ' ' )
    {
      if ( state == 4 )
      {
        state = 5;
      }
      else if ( state == 6 )
      {
        state = 7;
      }
      else if ( state == 10 )
      {
        state = 11;
      }
      else if ( state == 12 )
      {
        state = 13;
      }
    }
//
//  LEFT PAREN.
//
    else if ( c == '(' )
    {
      if ( state < 2 )
      {
        paren_sum = paren_sum + LEFT;
      }
      else
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current state = " << state << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }
    }
//
//  DIGIT (R, F, or W)
//
    else if ( ch_is_digit ( c ) )
    {
      if ( state <= 3 )
      {
        state = 4;
        *r = ch_to_digit ( c );
      }
      else if ( state == 4 )
      {
        d = ch_to_digit ( c );
        *r = 10 * (*r) + d;
      }
      else if ( state == 6 || state == 7 )
      {
        if ( *code == '*' )
        {
          if ( debug )
          {
            cerr << "\n";
            cerr << "S_TO_FORMAT - Fatal error!\n";
            cerr << "  Current state = " << state << "\n";
            cerr << "  Current code = '" << *code << "'.\n";
            cerr << "  Input character = '" << c << "'.\n";
            exit ( 1 );
          }
          state = -1;
          break;
        }
        state = 8;
        *w = ch_to_digit ( c );
      }
      else if ( state == 8 )
      {
        d = ch_to_digit ( c );
        *w = 10 * (*w) + d;
      }
      else if ( state == 9 )
      {
        state = 10;
        *m = ch_to_digit ( c );
      }
      else if ( state == 10 )
      {
        d = ch_to_digit ( c );
        *m = 10 * (*m) + d;
      }
      else
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current state = " << state << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }
    }
//
//  DECIMAL POINT.
//
    else if ( c == '.' )
    {
      if ( state == 8 )
      {
        state = 9;
      }
      else
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current state = " << state << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }
    }
//
//  RIGHT PAREN.
//
    else if ( c == ')' )
    {
      paren_sum = paren_sum + RIGHT;

      if ( paren_sum != 0 )
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current paren sum = " << paren_sum << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }

      if ( state == 6 && *code == '*' )
      {
        state = 12;
      }
      else if ( 6 <= state )
      {
        state = 12;
      }
      else
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current state = " << state << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }
    }
//
//  Format code
//
    else if ( ch_is_format_code ( c ) )
    {
      if ( state < 6 )
      {
        state = 6;
        *code = c;
      }
      else
      {
        if ( debug )
        {
          cerr << "\n";
          cerr << "S_TO_FORMAT - Fatal error!\n";
          cerr << "  Current state = " << state << "\n";
          cerr << "  Input character = '" << c << "'.\n";
          exit ( 1 );
        }
        state = -1;
        break;
      }
    }
//
//  Unexpected character
//
    else
    {
      if ( debug )
      {
        cerr << "\n";
        cerr << "S_TO_FORMAT - Fatal error!\n";
        cerr << "  Current state = " << state << "\n";
        cerr << "  Input character = '" << c << "'.\n";
        exit ( 1 );
      }
      state = -1;
      break;
    }
  }

  if ( paren_sum != 0 )
  {
    cerr << "\n";
    cerr << "S_TO_FORMAT - Fatal error!\n";
    cerr << "  Parentheses mismatch.\n";
    exit ( 1 );
  }

  if ( state < 0 )
  {
    cerr << "\n";
    cerr << "S_TO_FORMAT - Fatal error!\n";
    cerr << "  Parsing error.\n";
    exit ( 1 );
  }

  if ( *r == 0 )
  {
    *r = 1;
  }

  return;
}
//****************************************************************************80

int s_to_i4 ( string s, int &last, bool &error )

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
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int &LAST, the last character of S used to make IVAL.
//
//    Output, bool &ERROR is TRUE if an error occurred.
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

  last = 0;
  error = false;

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
        error = true;
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
        error = true;
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
        last = i - 1;
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
    last = s_len_trim ( s );
  }
  else
  {
    error = true;
    last = 0;
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
    ivec[i] = s_to_i4 ( s.substr(begin,length), lchar, error );

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

bool s_to_l ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_L reads an L from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Output, bool S_TO_L, the logical value.
//
{
  int i;
  bool l;
  int length;

  length = s.length ( );

  if ( length < 1 )
  {
    cerr << "\n";
    cerr << "S_TO_L - Fatal error!\n";
    cerr << "  Input string is empty.\n";
    exit ( 1 );
  }

  for ( i = 0; i < length; i++ )
  {
    if ( s[i] == '0' ||
         s[i] == 'f' ||
         s[i] == 'F' )
    {
      l = false;
      return l;
    }
    else if ( s[i] == '1' ||
              s[i] == 't' ||
              s[i] == 'T' )
    {
      l = true;
      return l;
    }
  }
  cerr << "\n";
  cerr << "S_TO_L - Fatal error!\n";
  cerr << "  Input did not contain boolean data.\n";
  exit ( 1 );
}
//****************************************************************************80

float s_to_r4 ( string s, int &lchar, bool &error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R4 reads an R4 from a string.
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
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2011
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
//    Output, int &LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool &ERROR, is true if an error occurred.
//
//    Output, float S_TO_R4, the real value that was read from the string.
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
  float r;
  float rbot;
  float rexp;
  float rtop;
  char TAB = 9;
  static float ten = 10.0;

  lchar = -1;
  error = false;

  nchar = s_len_trim ( s );
  r = 0.0;
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
    c = s[lchar+1];
    lchar = lchar + 1;
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
        lchar = lchar + 1;
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
        rtop = 10.0 * rtop + ( float ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( float ) ndig;
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
    if ( iterm == 1 || nchar <= lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && lchar + 1 == nchar )
  {
    lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    error = true;
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
      rexp = pow ( ten, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( ten, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r4vec ( string s, int n, float rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R4VEC reads an R4VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2011
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
//    Output, float RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R4VEC, is true if an error occurred.
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
    rvec[i] = s_to_r4 ( s.substr(begin,length), lchar, error );

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

double s_to_r8 ( string s, int &lchar, bool &error )

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
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2011
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
//    Output, int &LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool &ERROR, is true if an error occurred.
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
  static double ten = 10.0;

  lchar = -1;
  error = false;

  nchar = s_len_trim ( s );
  r = 0.0;
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
    c = s[lchar+1];
    lchar = lchar + 1;
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
        lchar = lchar + 1;
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
    if ( iterm == 1 || nchar <= lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && lchar + 1 == nchar )
  {
    lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    error = true;
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
      rexp = pow ( ten, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( ten, rexp );
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
    rvec[i] = s_to_r8 ( s.substr(begin,length), lchar, error );

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

string s_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_TRIM promotes the final null forward through trailing blanks.
//
//  Discussion:
//
//    What we're trying to say is that we reposition the null character
//    so that trailing blanks are no longer visible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be trimmed.
//
//    Output, string S_TRIM, the trimmed string.
//
{
  int i;
  int s_length;
  int s_nonblank_length;
  char *s2;
  string s3;

  s_length = s.length ( );

  s_nonblank_length = 0;
  for ( i = s_length - 1; 0 <= i; i-- )
  {
    if ( s[i] != ' ' )
    {
      s_nonblank_length = i + 1;
      break;
    }
  }

  s2 = new char[s_nonblank_length+1];

  for ( i = 0; i < s_nonblank_length; i++ )
  {
    s2[i] = s[i];
  }
  s2[s_length] = '\0';

  s3 = string ( s2 );

  delete [] s2;

  return s3;
}
//****************************************************************************80

string s_word_cap ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_CAP capitalizes the first character of each word in a string.
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
//    Input, string S, the string to be capitalized.
//
//    Output, string S_WORD_CAP, the word-capitalized string.
//
{
  bool blank;
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  blank = true;

  for ( i = 0; i < s_length; i++ )
  {
    if ( blank )
    {
      s2[i] = ch_cap ( s2[i] );
    }
    else
    {
      s2[i] = ch_low ( s2[i] );
    }
    blank = ( s2[i] == ' ' );
  }
  return s2;
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

void s_word_extract_first ( string s, string &s1, string &s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_EXTRACT_FIRST extracts the first word from a string.
//
//  Discussion:
//
//    A "word" is a string of characters terminated by a blank or
//    the end of the string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string.
//
//    Output, string &S1, the first word (initial blanks removed).
//
//    Output, string &S2, the remainder of the string, after removing
//    the first word (initial blanks removed).
//
{
  int i;
  int mode;
  int s_len;

  s_len = s.length ( );
  s1 = "";
  s2 = "";
  mode = 1;

  for ( i = 0; i < s_len; i++ )
  {
    if ( mode == 1 )
    {
      if ( s[i] != ' ' )
      {
         mode = 2;
      }
    }
    else if ( mode == 2 )
    {
      if ( s[i] == ' ' )
      {
        mode = 3;
      }
    }
    else if ( mode == 3 )
    {
      if ( s[i] != ' ' )
      {
        mode = 4;
      }
    }
    if ( mode == 2 )
    {
      s1 = s1 + s[i];
    }
    else if ( mode == 4 )
    {
      s2 = s2 + s[i];
    }
  }

  return;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ;; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

float swap_bytes_float ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    SWAP_BYTES_FLOAT swaps pairs of bytes in a float.
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
//    Input, float X, a float whose bytes are to be swapped.
//
//    Output, float SWAP_BYTES_FLOAT, a float with swapped bytes.
//
{
  unsigned char c;
  union
  {
    float yfloat;
    unsigned char ychar[4];
  } y;

  y.yfloat = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[1];
  y.ychar[1] = c;

  c = y.ychar[2];
  y.ychar[2] = y.ychar[3];
  y.ychar[3] = c;

  return ( y.yfloat );
}
//****************************************************************************80

int swap_bytes_int ( int x )

//****************************************************************************80
//
//  Purpose:
//
//    SWAP_BYTES_INT swaps pairs of bytes in an int.
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
//    Input, int X, an int whose bytes are to be swapped.
//
//    Output, int SWAP_BYTES_INT, an int with swapped bytes.
//
{
  char c;
  union
  {
    int yint;
    char ychar[4];
  } y;

  y.yint = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[1];
  y.ychar[1] = c;

  c = y.ychar[2];
  y.ychar[2] = y.ychar[3];
  y.ychar[3] = c;

  return ( y.yint );
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
//****************************************************************************80

string word_next_read ( string s, bool *done )

//****************************************************************************80
//
//  Purpose:
//
//    WORD_NEXT_READ "reads" words from a string, one at a time.
//
//  Discussion:
//
//    This routine was written to process tokens in a file.
//    A token is considered to be an alphanumeric string delimited
//    by whitespace, or any of various "brackets".
//
//    The following characters are considered to be a single word,
//    whether surrounded by spaces or not:
//
//      " ( ) { } [ ]
//
//    Also, if there is a trailing comma on the word, it is stripped off.
//    This is to facilitate the reading of lists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string, presumably containing words
//    separated by spaces.
//
//    Input/output, bool *DONE.
//    On input with a fresh string, set DONE to TRUE.
//    On output, the routine sets DONE:
//      FALSE if another word was read,
//      TRUE if no more words could be read.
//
//    Output, string WORD_NEXT_READ.
//    If DONE is FALSE, then WORD contains the "next" word read.
//    If DONE is TRUE, then WORD is NULL, because there was no more to read.
//
{
  int i;
  int ilo;
  int j;
  static int lenc = 0;
  static int next = 0;
  char TAB = 9;
  string word;
  char *word_chstar;
//
//  We "remember" LENC and NEXT from the previous call.
//
//  An input value of DONE = TRUE signals a new line of text to examine.
//
  if ( *done )
  {
    next = 0;
    *done = false;
    lenc = s.length ( );
    if ( lenc <= 0 )
    {
      *done = true;
      word = "\n";;
      return word;
    }
  }
//
//  Beginning at index NEXT, search the string for the next nonblank,
//  which signals the beginning of a word.
//
  ilo = next;
//
//  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
//
  for ( ; ; )
  {
    if ( lenc < ilo )
    {
      word = "\n";
      *done = true;
      next = lenc + 1;
      return word;
    }
//
//  If the current character is blank, skip to the next one.
//
    if ( s[ilo] != ' ' && s[ilo] != TAB )
    {
      break;
    }
    ilo = ilo + 1;
  }
//
//  ILO is the index of the next nonblank character in the string.
//
//  If this initial nonblank is a special character,
//  then that's the whole word as far as we're concerned,
//  so return immediately.
//
  if ( s[ilo] == '"' )
  {
    word = """";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '(' )
  {
    word = "(";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ')' )
  {
    word = ")";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '{' )
  {
    word = "{";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '}' )
  {
    word = "}";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '[' )
  {
    word = "[";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ']' )
  {
    word = "]";
    next = ilo + 1;
    return word;
  }
//
//  Now search for the last contiguous character that is not a
//  blank, TAB, or special character.
//
  next = ilo + 1;

  while ( next <= lenc )
  {
    if ( s[next] == ' ' )
    {
      break;
    }
    else if ( s[next] == TAB )
    {
      break;
    }
    else if ( s[next] == '"' )
    {
      break;
    }
    else if ( s[next] == '(' )
    {
      break;
    }
    else if ( s[next] == ')' )
    {
      break;
    }
    else if ( s[next] == '{' )
    {
      break;
    }
    else if ( s[next] == '}' )
    {
      break;
    }
    else if ( s[next] == '[' )
    {
      break;
    }
    else if ( s[next] == ']' )
    {
      break;
    }

    next = next + 1;
  }
//
//  Allocate WORD, copy characters, and return.
//
  if ( s[next-1] == ',' )
  {
    word_chstar = new char[next-ilo];
    i = 0;
    for ( j = ilo; j <= next - 2; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }
  else
  {
    word_chstar = new char[next+1-ilo];
    i = 0;
    for ( j = ilo; j <= next-1; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }

  return word;
}
