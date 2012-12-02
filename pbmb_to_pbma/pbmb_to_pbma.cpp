# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
void pbma_write ( string file_out_name, int xsize, int ysize, int *b );
void pbma_write_data ( ofstream &file_out, int xsize, int ysize, int *b );
void pbma_write_header ( ofstream &file_out, string file_out_name, int xsize,
  int ysize );
bool pbmb_check_data ( int xsize, int ysize, int *barray );
bool pbmb_read ( string input_name, int *xsize, int *ysize, int **barray );
bool pbmb_read_data ( ifstream &input, int xsize, int ysize, int *barray );
bool pbmb_read_header ( ifstream &input, int *xsize, int *ysize );
bool pbmb_to_pbma ( string input_name, string output_name );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PBMB_TO_PBMA.
//
//  Discussion:
//
//    PBMB_TO_PBMA converts a binary PBM file to an ASCII PBM file.
//
//    The application requires the PBMA_IO and PBMB_IO libraries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    pbmb_to_pbma file.pbmb file.pbma
//
//  Parameters:
//
//    FILE.PBMB is the name of the input binary PBM file to be read.
//
//    FILE.PBMA is the name of the output ASCII PBM file to be created.
//
{
  bool error;
  string input_name;
  string output_name;
  bool verbose = false;

  if ( verbose )
  {
    timestamp ( );
    cout << "\n";
    cout << "PBMB_TO_PBMA:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Convert a binary PBM file to ASCII PBM format.\n";
  }
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA:\n";
    cout << "  Please enter the name of the input binary PBM file:\n";

    cin >> input_name;

    error = cin.eof ( );

    if ( error )
    {
      exit ( 1 );
    }
  }
  else
  {
    input_name = argv[1];
  }
//
//  Get the specification for the output file.
//
  if ( argc < 3 )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA:\n";
    cout << "  Please enter the name of the output ASCII PBM file:\n";

    cin >> output_name;

    error = cin.eof( );

    if ( error )
    {
      exit ( 1 );
    }
  }
  else
  {
    output_name = argv[2];
  }

  error = pbmb_to_pbma ( input_name, output_name );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  if ( verbose )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA:\n";
    cout << "  Normal end of execution.\n";

    cout << "\n";
    timestamp ( );
  }

  return 0;
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

void pbma_write ( string file_out_name, int xsize, int ysize, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_WRITE writes the header and data for an ASCII PBM file.
//
//  Example:
//
//    P1
//    # feep.pbm
//    24 7
//    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
//    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
//    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
//    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
//    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
//    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_OUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *B, the array of XSIZE by YSIZE data values.
//
{
  ofstream file_out;
  int i;
  int *indexg;
  int j;
//
//  Open the output file.
//
  file_out.open ( file_out_name.c_str ( ) );

  if ( !file_out )
  {
    cerr << "\n";
    cerr << "PBMA_WRITE - Fatal error!\n";
    cerr << "  Cannot open the output file \"" << file_out_name << "\".\n";
    exit ( 1 );
  }
//
//  Write the header.
//
  pbma_write_header ( file_out, file_out_name, xsize, ysize );
//
//  Write the data.
//
  pbma_write_data ( file_out, xsize, ysize, b );
//
//  Close the file.
//
  file_out.close ( );

  return;
}
//****************************************************************************80

void pbma_write_data ( ofstream &file_out, int xsize, int ysize, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_WRITE_DATA writes the data for an ASCII PBM file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *B, the arrays of XSIZE by YSIZE data values.
//
{
  int i;
  int *indexb;
  int j;
  int numval;

  indexb = b;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      file_out << *indexb << " ";
      numval = numval + 1;
      indexb = indexb + 1;

      if ( numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize )
      {
        file_out << "\n";
      }
      else
      {
        file_out << " ";
      }

    }
  }
  return;
}
//****************************************************************************80

void pbma_write_header ( ofstream &file_out, string file_out_name, int xsize,
  int ysize )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_WRITE_HEADER writes the header of an ASCII PBM file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file.
//
//    Input, string FILE_OUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
{
  file_out << "P1\n";
  file_out << "# " << file_out_name << " created by PBMA_IO::PBMA_WRITE.\n";
  file_out << xsize << "  " << ysize << "\n";

  return;
}
//****************************************************************************80

bool pbmb_check_data ( int xsize, int ysize, int *barray )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_CHECK_DATA checks the data for a binary portable bit map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *BARRAY, the array of XSIZE by YSIZE bits.
//
//    Output, bool PBMB_CHECK_DATA, is true if an error occurred.
//
{
  int i;
  int *indexb;
  int j;

  if ( xsize <= 0 )
  {
    cout << "\n";
    cout << "PBMB_CHECK_DATA - Fatal error!\n";
    cout << "  XSIZE <= 0\n";
    cout << "  XSIZE = " << xsize << "\n";
    return true;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PBMB_CHECK_DATA - Fatal error!\n";
    cout << "  YSIZE <= 0\n";
    cout << "  YSIZE = " << ysize << "\n";
    return true;
  }

  if ( barray == NULL )
  {
    cout << "\n";
    cout << "PBMB_CHECK_DATA - Fatal error!\n";
    cout << "  Null pointer to data.\n";
    return true;
  }

  indexb = barray;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *indexb != 0 && *indexb != 1 )
      {
        cout << "\n";
        cout << "PBMB_CHECK_DATA - Fatal error!\n";
        cout << "  b(" << i << "," << j << ") = "
          << *indexb << ".\n";
        return true;
      }

      indexb = indexb + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool pbmb_read ( string input_name, int *xsize, int *ysize, int **barray )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_READ reads the header and data from a binary portable bit map file.
// 
//  Discussion:
//
//    Thanks to Jonas Schwertfeger for pointing out that, especially on Microsoft
//    Windows systems, a binary file needs to be opened as a binary file!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    02 April 2005
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_NAME, the name of the file containing the binary
//    portable bit map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int **BARRAY, the array of XSIZE by YSIZE data values.
//
//    Output, bool PBMB_READ, is true if an error occurred.
//
{
  bool error;
  ifstream input;

  input.open ( input_name.c_str ( ), ios::binary );

  if ( !input )
  {
    cout << "\n";
    cout << "PBMB_READ: Fatal error!\n";
    cout << "  Cannot open the input file " << input_name << "\n";
    return true;
  }
//
//  Read the header.
//
  error = pbmb_read_header ( input, xsize, ysize );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_READ: Fatal error!\n";
    cout << "  PBMB_READ_HEADER failed.\n";
    return true;
  }
//
//  Allocate storage for the data.
//
  *barray = new int [ (*xsize) * (*ysize) ];
//
//  Read the data.
//
  error = pbmb_read_data ( input, *xsize, *ysize, *barray );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_READ: Fatal error!\n";
    cout << "  PBMB_READ_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  input.close ( );

  return false;
}
//****************************************************************************80

bool pbmb_read_data ( ifstream &input, int xsize, int ysize, int *barray )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_READ_DATA reads the data in a binary portable bit map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the file containing the binary
//    portable bit map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *BARRAY, the array of XSIZE by YSIZE data values.
//
//    Output, bool PBMB_READ_DATA, is true if an error occurred.
//
{
  int bit;
  char c;
  unsigned char c2;
  int i;
  int *indexb;
  int j;
  int k;
  int numbyte;

  indexb = barray;
  numbyte = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( i%8 == 0 )
      {
        input.read ( &c, 1 );
        c2 = ( unsigned char ) c;
        if ( input.eof ( ) )
        {
          cout << "\n";
          cout << "PBMB_CHECK_DATA - Fatal error!\n";
          cout << "  Failed reading byte " << numbyte << "\n";
          return true;
        }
        numbyte = numbyte + 1;
      }

      k = 7 - i%8;
      bit = ( c2 >> k )%2;

      *indexb = bit;
      indexb = indexb + 1;
    }
  }
  return false;
}
//****************************************************************************80

bool pbmb_read_header ( ifstream &input, int *xsize, int *ysize )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_READ_HEADER reads the header of a binary portable bit map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the file containing the binary
//    portable bit map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, bool PBMB_READ_HEADER, is true if an error occurred.
//
{
  int count;
  string line;
  string rest;
  int step;
  int width;
  string word;

  step = 0;

  while ( 1 )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "PBMB_READ_HEADER - Fatal error!\n";
      cout << "  End of file.\n";
      return true;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    if ( step == 0 )
    {
      s_word_extract_first ( line, word, rest );

      if ( s_len_trim ( word ) <= 0 )
      {
        continue;
      }
      
      if ( !s_eqi ( word, "P4" ) )
      {
        cout << "\n";
        cout << "PBMB_READ_HEADER - Fatal error.\n";
        cout << "  Bad magic number = \"" << word << "\".\n";
        return true;
      }
      line = rest;
      step = 1;
    }

    if ( step == 1 )
    {
      s_word_extract_first ( line, word, rest );
 
      if ( s_len_trim ( word ) <= 0 )
      {
        continue;
      }
      *xsize = atoi ( word.c_str ( ) );
      line = rest;
      step = 2;
    }

    if ( step == 2 )
    {
      s_word_extract_first ( line, word, rest );
 
      if ( s_len_trim ( word ) <= 0 )
      {
        continue;
      }
      *ysize = atoi ( word.c_str ( ) );
      break;
    }

  }

  return false;
}
//****************************************************************************80

bool pbmb_to_pbma ( string input_name, string output_name )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_TO_PBMA converts one PBMB file to PBMA format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_NAME, the name of the input binary PBM file.
//
//    Input, string OUTPUT_NAME, the name of the output ASCII PBM file.
//
//    Output, bool HANDLE, is true if an error occurred.
//
{
  int *b;
  int *b2;
  bool error;
  int xsize;
  int ysize;
//
//  Read the input file.
//
  error = pbmb_read ( input_name, &xsize, &ysize, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA: Fatal error!\n";
    cout << "  PBMB_READ failed.\n";
    return true;
  }
//
//  Check the data.
//
  error = pbmb_check_data ( xsize, ysize, b );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_TO_PBMA: Fatal error!\n";
    cout << "  PBMB_CHECK_DATA reports bad data from the file.\n";

    delete [] b;
    return true;
  }
//
//  Write the output file.
//
  pbma_write ( output_name, xsize, ysize, b );

  delete [] b;

  return false;
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
