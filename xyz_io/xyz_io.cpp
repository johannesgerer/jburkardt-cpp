# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <fstream>

using namespace std;

# include "xyz_io.hpp"

//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

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
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 ) 
  {
    c1 = c1 - 32;
  } 
  if ( 97 <= c2 && c2 <= 122 ) 
  {
    c2 = c2 - 32;
  }     

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
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
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
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

void i4vec_copy ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_COPY copies an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, int A1[N], the vector to be copied.
//
//    Output, int A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
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

void timestamp ( void )

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
//    24 September 2003
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
//****************************************************************************80

void xyz_data_print ( int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_DATA_PRINT prints the data for an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double XY[3*POINT_NUM], the arrays of coordinate data.
//
{
  int j;

  cout << "\n";
  for ( j = 0; j < point_num; j++ )
  {
    cout << setw(10) << xyz[0+j*3] << "  "
         << setw(10) << xyz[1+j*3] << "  " 
         << setw(10) << xyz[2+j*3] << "\n";
    
  }
  return;
}
//****************************************************************************80

void xyz_data_read ( string input_filename, int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_DATA_READ reads the data in an XYZ file.
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
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double XYZ[3*POINT_NUM], the point coordinates.
//
{
  bool error;
  int i;
  ifstream input;
  int j;
  string text;
  double temp[3];

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZ_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  j = 0;

  while ( j < point_num )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }
    
    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }
//
//  Extract two real numbers.
//
    error = s_to_r8vec ( text, 3, temp );

    if ( error )
    {
      cout << "\n";
      cout << "XYZ_DATA_READ - Fatal error!\n";
      cout << "  S_TO_R8VEC returned error flag.\n";
      exit ( 1 );
    }

    xyz[0+j*3] = temp[0];
    xyz[1+j*3] = temp[1];
    xyz[2+j*3] = temp[2];
    j = j + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyz_data_write ( ofstream &output_unit, int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_DATA_WRITE writes the data for an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double XY[3*POINT_NUM], the arrays of coordinate data.
//
{
  int j;

  for ( j = 0; j < point_num; j++ )
  {
    output_unit << setw(10) << xyz[0+j*3] << "  "
                << setw(10) << xyz[1+j*3] << "  " 
                << setw(10) << xyz[2+j*3] << "\n";
    
  }
  return;
}
//****************************************************************************80

void xyz_example ( int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_EXAMPLE sets up data suitable for an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Output, double XYZ[3*POINT_NUM], the point coordinates.
//
{
  double a;
  int i;
  double r;
  double r8_pi = 3.141592653589793;
  double theta;
  double theta1;
  double theta2;

  r = 1.0;
  theta1 = 0.0;
  theta2 = 3.0 * 2.0 * r8_pi;
  a = 2.0 / ( theta2 - theta1 );

  for ( i = 0; i < point_num; i++ )
  {
    if ( point_num == 1 )
    {
      theta = 0.5 * ( theta1 + theta2 );
    }
    else
    {
      theta = ( ( double ) ( point_num - i - 1 ) * theta1 
              + ( double ) (             i     ) * theta2 ) 
              / ( double ) ( point_num     - 1 );
    }

    xyz[0+i*3] = r * cos ( theta );
    xyz[1+i*3] = r * sin ( theta );
    xyz[2+i*3] = a * theta;
  }

  return;
}
//****************************************************************************80

int xyz_example_size ( )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_EXAMPLE_SIZE sizes an example XYZ dataset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, integer XYZ_EXAMPLE_SIZE, the number of points.
//
{
  int point_num;

  point_num = 101;

  return point_num;
}
//****************************************************************************80

void xyz_header_print ( int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_HEADER_PRINT prints the header of an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
{
  cout << "\n";
  cout << "  Number of points = " << point_num << "\n";

  return;
}
//****************************************************************************80

void xyz_header_read ( string input_filename, int *point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_HEADER_READ reads the header of an XYZ file.
//
//  Discussion:
//
//    All we do here is count the number of lines that are not comments 
//    and not blank.  Each such line is assumed to represent a single point 
//    coordinate record.
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
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Output, int *POINT_NUM, the number of points.
//
{
  ifstream input;
  string text;

  *point_num = 0;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZ_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  while ( 1 )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    *point_num = *point_num + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyz_header_write ( string output_filename, ofstream &output_unit, 
  int point_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_HEADER_WRITE writes the header of an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
//
//    Input, int POINT_NUM, the number of points.
//
{
  output_unit << "#  " << output_filename << "\n";
  output_unit << "#  created by xyz_io::xyz_header_write.C\n";
  output_unit << "#\n";
  output_unit << "#  Number of points = " << point_num << "\n";
  output_unit << "#\n";

  return;
}
//****************************************************************************80

void xyz_read ( string input_filename, int *point_num, double *xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_READ reads the header and data from an XYZ file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    05 January 2009
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Output, int *POINT_NUM, the number of points.
//
//    Output, double *XYZ[3*(*POINT_NUM)], the point coordinates.
//
{
  xyz_header_read ( input_filename, point_num );

  *xyz = new double[3 * (*point_num) ];

  xyz_data_read ( input_filename, *point_num, *xyz );

  return;
}
//****************************************************************************80

void xyz_read_test ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_READ_TEST tests the XYZ file read routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
{
  int point_num;
  double *xyz = NULL;

  xyz_read ( input_filename, &point_num, &xyz );

  delete [] xyz;

  return;
}
//****************************************************************************80

void xyz_write ( string output_filename, int point_num, double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_WRITE writes the header and data for an XYZ file.
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    01 January 2009
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double XYZ[3*POINT_NUM], the arrays of coordinate data.
//
{
  ofstream output_unit;
//
//  Open the output file.
//
  output_unit.open ( output_filename.c_str() );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "XYZ_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }
//
//  Write the header.
//
  xyz_header_write ( output_filename, output_unit, point_num );
//
//  Write the data.
//
  xyz_data_write ( output_unit, point_num, xyz );
//
//  Close the file.
//
  output_unit.close ( );

  return;
}
//****************************************************************************80

void xyz_write_test ( string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    XYZ_WRITE_TEST tests the XYZ write routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file to contain the data.
//
{
  int point_num;
  double *xyz;
//
//  Set the data.
//
  point_num = xyz_example_size ( );
//
//  Allocate memory.
//
  xyz = new double[3*point_num];
//
//  Set the data.
//
  xyz_example ( point_num, xyz );
//
//  Write the data to the file.
//
  xyz_write ( output_filename, point_num, xyz );

  delete [] xyz;
 
  return;
}
//****************************************************************************80

void xyzf_data_print ( int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_DATA_PRINT prints the data of an XYZF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Input, int FACE_DATA[FACE_DATA_NUM], indices
//    of points that form faces.
//
{
  int i;
  int face;

  cout << "\n";
  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(4) << face << "  "
         << "  " << setw(8) << face_pointer[face]
         << "  " << setw(8) << face_pointer[face+1] - 1 << "\n"; 
  }

  cout << "\n";
  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(4) << face << "  ";
    for ( i = face_pointer[face]; i < face_pointer[face+1]; i++ )
    {
      cout << "  " << setw(4) << face_data[i];
    }
    cout << "\n";    
  }
  return;
}
//****************************************************************************80

void xyzf_data_read ( string input_filename, int face_num, int face_data_num,
  int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_DATA_READ reads the data in an XYZF file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    FACE ITEMS, which are indices of points on a face.
//
//    The routine ignores comments and blanks and returns
//    the number of face items.
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
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Output, int FACE_DATA[FACE_DATA_NUM], the face items.
//
{
  int face;
  int ierror;
  int ilo;
  ifstream input;
  int n;
  string text;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZF_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  face = 0;
  face_pointer[0] = 0;

  while ( face < face_num )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "XYZF_DATA_READ - Fatal error!\n";
      cout << "  Unexpected end of information.\n";
      exit ( 1 );
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );
    face_pointer[face+1] = face_pointer[face] + n;
 
    ilo = face_pointer[face];

    ierror = s_to_i4vec ( text, n, face_data+ilo );

    if ( ierror != 0 )
    {
      cout << "\n";
      cout << "XYZF_DATA_READ - Fatal error!\n";
      cout << "  Error from S_TO_I4VEC.\n";
      exit ( 1 );
    }
    face = face + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyzf_data_write ( ofstream &output_unit, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_DATA_WRITE writes the data of an XYZF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the XY file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Input, int FACE_DATA[FACE_DATA_NUM], indices
//    of points that form faces.
//
{
  int i;
  int face;

  for ( face = 0; face < face_num; face++ )
  {
    for ( i = face_pointer[face]; i < face_pointer[face+1]; i++ )
    {
      output_unit << "  " << face_data[i];
    }
    output_unit << "\n";    
  }
  return;
}
//****************************************************************************80

void xyzf_example ( int point_num, int face_num, int face_data_num, 
  double xyz[], int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_EXAMPLE sets data suitable for a pair of XYZ and XYZF files.
//
//  Discussion:
//
//    There are 8 points.
//    There are 6 faces.
//    There are 24 face items.
//
//       8------7
//      /|     /|
//     / |    / |
//    5------6  |
//    |  4---|--3
//    | /    | /
//    |/     |/
//    1------2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Output, double XYZ[3*POINT_NUM], the point coordinates.
//
//    Output, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Output, int FACE_DATA[FACE_DATA_NUM], indices
//    of points that form faces.
//
{
# define FACE_DATA_NUM 24
# define FACE_NUM 6
# define POINT_NUM 8

  int FACE_DATA[FACE_DATA_NUM] = {
     0, 3, 2, 1, 
     1, 2, 6, 5, 
     4, 5, 6, 7, 
     4, 7, 3, 0, 
     0, 1, 5, 4, 
     2, 3, 7, 6 };
  int FACE_POINTER[FACE_NUM+1] = { 0, 4, 8, 12, 16, 20, 24 };
  double XYZ[3*POINT_NUM] = {
     0.0,  0.0,  0.0, 
     1.0,  0.0,  0.0, 
     1.0,  1.0,  0.0, 
     0.0,  1.0,  0.0, 
     0.0,  0.0,  1.0, 
     1.0,  0.0,  1.0, 
     1.0,  1.0,  1.0, 
     0.0,  1.0,  1.0 };

  i4vec_copy ( face_data_num, FACE_DATA, face_data );
  i4vec_copy ( face_num + 1, FACE_POINTER, face_pointer );
  r8vec_copy ( 3 * point_num, XYZ, xyz );

  return;
# undef FACE_DATA_NUM
# undef FACE_NUM
# undef POINT_NUM
}
//****************************************************************************80

void xyzf_example_size ( int *point_num, int *face_num, int *face_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_EXAMPLE_SIZE sizes the data to be created by XYZF_EXAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *POINT_NUM, the number of points.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_DATA_NUM, the number of face items.
//
{
  *face_data_num = 24;
  *face_num = 6;
  *point_num = 8;

  return;
}
//****************************************************************************80

void xyzf_header_print ( int point_num, int face_num, int face_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_HEADER_PRINT prints the header of an XYZF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
{
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of faces      = " << face_num << "\n";
  cout << "  Number of face items = " << face_data_num << "\n";

  return;
}
//****************************************************************************80

void xyzf_header_read ( string input_filename, int *face_num, 
  int *face_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_HEADER_READ determines the number of face items in an XYZF file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    FACE ITEMS, which are indices of points on a face.
//
//    The routine ignores comments and blanks and returns
//    the number of face items.
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
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *FACE_NUM, the number of faces.
//
//    Output, int *FACE_DATA_NUM, the number of face items.
//
{
  int i;
  int i4_val;
  int ierror;
  ifstream  input;
  int length;
  int n;
  string text;

  *face_data_num = 0;
  *face_num = 0;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZF_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );

    *face_data_num = *face_data_num + n;

    *face_num = *face_num + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyzf_header_write ( string output_filename, ofstream &output_unit, 
  int point_num, int face_num, int face_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_HEADER_WRITE writes the header of an XYZF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
{
  output_unit << "#  " << output_filename << "\n";
  output_unit << "#  created by xyz_io::xyzf_header_write.C\n";
  output_unit << "#\n";
  output_unit << "#  Number of points     = " << point_num << "\n";
  output_unit << "#  Number of faces      = " << face_num << "\n";
  output_unit << "#  Number of face items = " << face_data_num << "\n";
  output_unit << "#\n";

  return;
}
//****************************************************************************80

void xyzf_write ( string output_filename, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZF_WRITE writes the header and data for an XYZF file.
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    07 January 2009
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_DATA_NUM, the number of face items.
//
//    Input, int FACE_POINTER[FACE_NUM+1], pointers to the
//    first face item for each face.
//
//    Input, int FACE_DATA[FACE_DATA_NUM], indices
//    of points that form faces.
//
{
  ofstream output_unit;
//
//  Open the output file.
//
  output_unit.open ( output_filename.c_str() );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "XYZF_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }
//
//  Write the header.
//
  xyzf_header_write ( output_filename, output_unit, point_num, face_num, 
    face_data_num );
//
//  Write the data.
//
  xyzf_data_write ( output_unit, point_num, face_num, face_data_num, face_pointer, 
    face_data );
//
//  Close the file.
//
  output_unit.close ( );

  return;
}
//****************************************************************************80

void xyzl_data_print ( int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_DATA_PRINT prints the data of an XYZL file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
//    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
//    first line item for each line.
//
//    Input, int LINE_DATA[LINE_DATA_NUM], indices
//    of points that form lines.
//
{
  int i;
  int line;

  cout << "\n";
  for ( line = 0; line < line_num; line++ )
  {
    cout << "  " << setw(4) << line << "  "
         << "  " << setw(8) << line_pointer[line]
         << "  " << setw(8) << line_pointer[line+1] - 1 << "\n"; 
  }

  cout << "\n";
  for ( line = 0; line < line_num; line++ )
  {
    cout << "  " << setw(4) << line << "  ";
    for ( i = line_pointer[line]; i < line_pointer[line+1]; i++ )
    {
      cout << "  " << setw(4) << line_data[i];
    }
    cout << "\n";    
  }
  return;
}
//****************************************************************************80

void xyzl_data_read ( string input_filename, int line_num, int line_data_num,
  int line_pointer[], int line_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_DATA_READ reads the data in an XYZL file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    LINE ITEMS, which are indices of points on a line.
//
//    The routine ignores comments and blanks and returns
//    the number of line items.
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
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
//    Output, int LINE_POINTER[LINE_NUM+1], pointers to the
//    first line item for each line.
//
//    Output, int LINE_DATA[LINE_DATA_NUM], the line items.
//
{
  int ierror;
  int ilo;
  ifstream input;
  int line;
  int n;
  string text;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZL_DATA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  line = 0;
  line_pointer[0] = 0;

  while ( line < line_num )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "XYZL_DATA_READ - Fatal error!\n";
      cout << "  Unexpected end of information.\n";
      exit ( 1 );
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );
    line_pointer[line+1] = line_pointer[line] + n;
 
    ilo = line_pointer[line];

    ierror = s_to_i4vec ( text, n, line_data+ilo );

    if ( ierror != 0 )
    {
      cout << "\n";
      cout << "XYZL_DATA_READ - Fatal error!\n";
      cout << "  Error from S_TO_I4VEC.\n";
      exit ( 1 );
    }
    line = line + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyzl_data_write ( ofstream &output_unit, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_DATA_WRITE writes the data of an XYZL file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the XY file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
//    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
//    first line item for each line.
//
//    Input, int LINE_DATA[LINE_DATA_NUM], indices
//    of points that form lines.
//
{
  int i;
  int line;

  for ( line = 0; line < line_num; line++ )
  {
    for ( i = line_pointer[line]; i < line_pointer[line+1]; i++ )
    {
      output_unit << "  " << line_data[i];
    }
    output_unit << "\n";    
  }
  return;
}
//****************************************************************************80

void xyzl_example ( int point_num, int line_num, int line_data_num, 
  double xyz[], int line_pointer[], int line_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_EXAMPLE sets data suitable for a pair of XYZ and XYZL files.
//
//  Discussion:
//
//    There are 8 points.
//    There are 6 lines.
//    There are 18 line items.
//
//       8------7
//      /|     /|
//     / |    / |
//    5------6  |
//    |  4---|--3
//    | /    | /
//    |/     |/
//    1------2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
//    Output, double XYZ[3*POINT_NUM], the point coordinates.
//
//    Output, int LINE_POINTER[LINE_NUM+1], pointers to the
//    first line item for each line.
//
//    Output, int LINE_DATA[LINE_DATA_NUM], indices
//    of points that form lines.
//
{
# define LINE_DATA_NUM 18
# define LINE_NUM 6
# define POINT_NUM 8

  int LINE_DATA[LINE_DATA_NUM] = {
     0,  1,  2,  3,  0, 
     4,  5,  6,  7,  4, 
     0,  4, 
     1,  5, 
     2,  6, 
     3,  7 };
  int LINE_POINTER[LINE_NUM+1] = { 0, 5, 10, 12, 14, 16, 18 };
  double XYZ[3*POINT_NUM] = {
     0.0,  0.0,  0.0, 
     1.0,  0.0,  0.0, 
     1.0,  1.0,  0.0, 
     0.0,  1.0,  0.0, 
     0.0,  0.0,  1.0, 
     1.0,  0.0,  1.0, 
     1.0,  1.0,  1.0, 
     0.0,  1.0,  1.0 };

  i4vec_copy ( line_data_num, LINE_DATA, line_data );
  i4vec_copy ( line_num + 1, LINE_POINTER, line_pointer );
  r8vec_copy ( 3 * point_num, XYZ, xyz );

  return;
# undef LINE_DATA_NUM
# undef LINE_NUM
# undef POINT_NUM
}
//****************************************************************************80

void xyzl_example_size ( int *point_num, int *line_num, int *line_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_EXAMPLE_SIZE sizes the data to be created by XYZL_EXAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *POINT_NUM, the number of points.
//
//    Output, int *LINE_NUM, the number of lines.
//
//    Output, int *LINE_DATA_NUM, the number of line items.
//
{
  *line_data_num = 18;
  *line_num = 6;
  *point_num = 8;

  return;
}
//****************************************************************************80

void xyzl_header_print ( int point_num, int line_num, int line_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_HEADER_PRINT prints the header of an XYZL file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
{
  cout << "\n";
  cout << "  Number of points     = " << point_num << "\n";
  cout << "  Number of lines      = " << line_num << "\n";
  cout << "  Number of line items = " << line_data_num << "\n";

  return;
}
//****************************************************************************80

void xyzl_header_read ( string input_filename, int *line_num, 
  int *line_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_HEADER_READ determines the number of line items in an XYZL file.
//
//  Discussion:
//
//    This routine assumes that the file contains exactly three kinds of
//    records:
//
//    COMMENTS which begin with a '#' character in column 1;
//    BLANKS which contain nothing but 'whitespace';
//    LINE ITEMS, which are indices of points on a line.
//
//    The routine ignores comments and blanks and returns
//    the number of line items.
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
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *LINE_NUM, the number of lines.
//
//    Output, int *LINE_DATA_NUM, the number of line items.
//
{
  int i;
  int i4_val;
  int ierror;
  ifstream  input;
  int length;
  int n;
  string text;

  *line_data_num = 0;
  *line_num = 0;

  input.open ( input_filename.c_str() );

  if ( !input )
  {
    cout << "\n";
    cout << "XYZL_HEADER_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_filename << "\".\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( text[0] == '#' || s_len_trim ( text ) == 0 )
    {
      continue;
    }

    n = s_word_count ( text );

    *line_data_num = *line_data_num + n;

    *line_num = *line_num + 1;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void xyzl_header_write ( string output_filename, ofstream &output_unit, 
  int point_num, int line_num, int line_data_num )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_HEADER_WRITE writes the header of an XYZL file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, ofstream &OUTPUT_UNIT, a pointer to the file to contain the data.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
{
  output_unit << "#  " << output_filename << "\n";
  output_unit << "#  created by xyz_io::xyzl_header_write.C\n";
  output_unit << "#\n";
  output_unit << "#  Number of points     = " << point_num << "\n";
  output_unit << "#  Number of lines      = " << line_num << "\n";
  output_unit << "#  Number of line items = " << line_data_num << "\n";
  output_unit << "#\n";

  return;
}
//****************************************************************************80

void xyzl_write ( string output_filename, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    XYZL_WRITE writes the header and data for an XYZL file.
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    05 January 2009
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the name of the file.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int LINE_NUM, the number of lines.
//
//    Input, int LINE_DATA_NUM, the number of line items.
//
//    Input, int LINE_POINTER[LINE_NUM+1], pointers to the
//    first line item for each line.
//
//    Input, int LINE_DATA[LINE_DATA_NUM], indices
//    of points that form lines.
//
{
  ofstream output_unit;
//
//  Open the output file.
//
  output_unit.open ( output_filename.c_str() );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "XYZL_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_filename << "\".\n";
    exit ( 1 );
  }
//
//  Write the header.
//
  xyzl_header_write ( output_filename, output_unit, point_num, line_num, 
    line_data_num );
//
//  Write the data.
//
  xyzl_data_write ( output_unit, point_num, line_num, line_data_num, line_pointer, 
    line_data );
//
//  Close the file.
//
  output_unit.close ( );

  return;
}
