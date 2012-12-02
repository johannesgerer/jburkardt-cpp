# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

# include "vector_read.hpp"

//****************************************************************************80

int i4vec_read ( ifstream &input, int **i4vec )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_READ reads the "next" I4VEC from the input stream.
//
//  Discussion:
//
//    Before calling this routine, you must open the input stream.
//    This means your file must include the header
//
//      # include <fstream>
//
//    Then your program would declare the input stream like this:
//
//      ifstream input;
//
//    and open an input file to be associated with the stream like this:
//
//      input.open ( input_filename );
//
//    and then pass a reference to the input stream like this:
//
//      n = i4vec_read ( &input, &i4vec );
//
//    If N is positive, then I4VEC is a pointer to an integer vector of N
//    data items.  Once the user is done with the data, it can be freed
//    by a command like:
//
//      delete [] i4vec;
//
//    Repeated calls to I4VEC_READ will result in successive lines of the
//    input file being read.  Once you are done, close the file by
//
//      input.close ( );
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the input stream.
//
//    Output, int **I4VEC, the address of the vector.
//
//    Output, int I4VEC_READ, the number of entries in the vector, or
//    -1 if no data could be read.
//
{
  bool error;
  char line[255];
  int n;

  n = 0;
//
//  First, read repeatedly until you get a line that is not blank
//  and which does not begin with a comment character.
//
  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      n = -1;
      *i4vec = NULL;
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) <= 0 )
    {
      continue;
    }
//
//  We have a line.  Count the number of "words" in it.
//
    n = s_word_count ( line );

    if ( n <= 0 )
    {
      continue;
    }

    *i4vec = new int[n];
//
//  Extract the values.
//
    error = s_to_i4vec ( line, n, *i4vec );

    if ( !error )
    {
      break;
    }
    delete [] *i4vec;
    n = 0;
  }
  return n;
}
//****************************************************************************80

int s_len_trim ( char *s )

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
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
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

int s_to_i4 ( char *s, int *last, bool *error )

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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a string to be examined.
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

  while ( *s )
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

bool s_to_i4vec ( char *s, int n, int i4vec[] )

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
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int I4VEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  bool error;
  int i;
  int lchar;

  error = false;

  for ( i = 0; i < n; i++ )
  {
    i4vec[i] = s_to_i4 ( s, &lchar, &error );

    if ( error )
    {
      cout << "\n";
      cout << "S_TO_I4VEC - Fatal error!\n";
      cout << "  S_TO_I4 returned error while reading item " << i << "\n";
      return error;
    }
    s = s + lchar;
  }

  return error;
}
//****************************************************************************80

int s_word_count ( char *s )

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
//    30 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int i;
  int word_num;

  word_num = 0;
  blank = true;

  while ( *s )
  {
    if ( *s == ' ' )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_num = word_num + 1;
      blank = false;
    }
    *s++;
  }

  return word_num;
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
