# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
void pgma_write ( string output_name, int xsize, int ysize, int *g );
void pgma_write_data ( ofstream &output, int xsize, int ysize, int *g );
void pgma_write_header ( ofstream &output, string output_name, int xsize, 
  int ysize, int maxg );
bool pgmb_check_data ( int xsize, int ysize, unsigned char maxg, 
  unsigned char *g );
bool pgmb_read ( string input_name, int *xsize, int *ysize, 
  unsigned char *maxg, unsigned char **g );
bool pgmb_read_data ( ifstream &input, int xsize, int ysize, 
  unsigned char *g );
bool pgmb_read_header ( ifstream &input, int *xsize, int *ysize, 
  unsigned char *maxg );
bool pgmb_to_pgma ( char *file_in_name, char *file_out_name );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( );
void ucvec_to_i4vec ( int n, unsigned char *a, int *b );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PGMB_TO_PGMA.
//
//  Discussion:
//
//    PGMB_TO_PGMA converts a binary PGM file to an ASCII PGM file.
//
//  Usage:
//
//    pgmb_to_pgma file.pgmb file.pgma
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
//  Parameters:
//
//    FILE.PGMB is the name of the input binary PGM file to be read.
//
//    FILE.PGMA is the name of the output ASCII PGM file to be created.
//  
{  
  bool error;
  char file_in_name[80];
  char file_out_name[80];
  bool verbose = false;

  if ( verbose )
  {
    timestamp ( );
    cout << "\n";
    cout << "PGMB_TO_PGMA:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Convert a binary PGM file to ASCII PGM format.\n";
  }
//
//  Get the specification for the input file.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA:\n";
    cout << "  Please enter the input PGMB file name:\n";
    
    cin.getline ( file_in_name, sizeof ( file_in_name ) );

    error = cin.eof();

    if ( error )
    {
      exit ( 1 );
    }
  }
  else 
  {
    strcpy ( file_in_name, argv[1] );
  }
//
//  Get the specification for the output file.
//
  if ( argc < 3 ) 
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA:\n";
    cout << "  Please enter the output PGMA file name:\n";

    cin.getline ( file_out_name, sizeof ( file_out_name ) );

    error = cin.eof();

    if ( error )
    {
      exit ( 1 );
    }
  }
  else 
  {
    strcpy ( file_out_name, argv[2] );
  }

  error = pgmb_to_pgma ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  if ( verbose )
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA:\n";
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

void pgma_write ( string output_name, int xsize, int ysize, int *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_WRITE writes the header and data for an ASCII PGM file.
// 
//  Example:
//
//    P2
//    # feep.pgm
//    24 7
//    15
//    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
//    0  3  3  3  3  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15 15 15 15  0
//    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0 15  0
//    0  3  3  3  0  0  0  7  7  7  0  0  0 11 11 11  0  0  0 15 15 15 15  0
//    0  3  0  0  0  0  0  7  0  0  0  0  0 11  0  0  0  0  0 15  0  0  0  0
//    0  3  0  0  0  0  0  7  7  7  7  0  0 11 11 11 11  0  0 15  0  0  0  0
//    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
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
//    Input, string OUTPUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *G, the array of XSIZE by YSIZE data values.
//
{
  ofstream output;
  int i;
  int *indexg;
  int j;
  int maxg;
//
//  Open the output file.
//
  output.open ( output_name.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "PGMA_WRITE - Fatal error!\n";
    cerr << "  Cannot open the output file \"" << output_name << "\".\n";
    exit ( 1 );
  }
//
//  Compute the maximum.
//
  maxg = 0;
  indexg = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( maxg < *indexg )
      {
        maxg = *indexg;
      }
      indexg = indexg + 1;

    }
  }
//
//  Write the header.
//
  pgma_write_header ( output, output_name, xsize, ysize, maxg );
//
//  Write the data.
//
  pgma_write_data ( output, xsize, ysize, g );
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void pgma_write_data ( ofstream &output, int xsize, int ysize, int *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_WRITE_DATA writes the data for an ASCII PGM file.
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
//    Input, ofstream &OUTPUT, a pointer to the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *G, the array of XSIZE by YSIZE data.
//
{
  int i;
  int *indexg;
  int j;
  int numval;

  indexg = g;
  numval = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      output << *indexg;
      numval = numval + 1;
      indexg = indexg + 1;

      if ( numval % 12 == 0 || i == xsize - 1 || numval == xsize * ysize )
      {
        output << "\n";
      }
      else
      {
        output << " ";
      }

    }
  }
  return;
}
//****************************************************************************80

void pgma_write_header ( ofstream &output, string output_name, int xsize, 
  int ysize, int maxg )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_WRITE_HEADER writes the header of an ASCII PGM file.
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
//    Input, ofstream &OUTPUT, a pointer to the file.
//
//    Input, string OUTPUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXG, the maximum gray value.
//
{
  output << "P2\n";
  output << "# " << output_name << " created by PGMA_IO::PGMA_WRITE.\n";
  output << xsize << "  " << ysize << "\n";
  output << maxg << "\n";

  return;
}
//****************************************************************************80

bool pgmb_check_data ( int xsize, int ysize, unsigned char maxg, 
  unsigned char *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_CHECK_DATA checks the data for a binary portable gray map file.
//
//  Discussion:
//
//    XSIZE and YSIZE must be positive, the pointers must not be null,
//    and the data must be nonnegative and no greater than MAXG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char MAXG, the maximum gray value.
//
//    Input, unsigned char *G, the array of XSIZE by YSIZE data values.
//
//    Output, bool PGMB_CHECK_DATA, is
//    true, if an error was detected, or
//    false, if the data was legal.
//
{
  int i;
  unsigned char *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cout << "\n";
    cout << "PGMB_CHECK_DATA: Error!\n";
    cout << "  xsize <= 0.\n";
    cout << "  xsize = " << xsize << "\n";
    return true;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PGMB_CHECK_DATA: Error!\n";
    cout << "  ysize <= 0.\n";
    cout << "  ysize = " << ysize << "\n";
    return true;
  }

  if ( g == NULL )
  {
    cout << "\n";
    cout << "PGMB_CHECK_DATA: Error!\n";
    cout << "  Null pointer to g.\n";
    return true;
  }

  index = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( maxg < *index )
      {
        cout << "\n";
        cout << "PGMB_CHECK_DATA - Fatal error!\n";
        cout << "  Data exceeds MAXG = " << ( int ) maxg << "\n";
        cout << "  G(" << i << "," << j << ")=" << ( int ) (*index) << "\n";
        return true;
      }

      index = index + 1;
    }
  } 

  return false;
}
//****************************************************************************80

bool pgmb_read ( string input_name, int *xsize, int *ysize, 
  unsigned char *maxg, unsigned char **g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_READ reads the header and data from a binary portable gray map file.
// 
//  Discussion:
//
//    Thanks to Jonas Schwertfeger for pointing out that, especially on 
//    Microsoft Windows systems, a binary file needs to be opened as a 
//    binary file!
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
//    portable gray map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, unsigned char *MAXG, the maximum gray value.
//
//    Output, unsigned char **G, the array of XSIZE by YSIZE data values.
//
//    Output, bool PGMB_READ, is true if an error occurred.
//
{
  bool error;
  ifstream input;
  int numbytes;

  input.open ( input_name.c_str ( ), ios::binary );

  if ( !input )
  {
    cout << "\n";
    cout << "PGMB_READ: Fatal error!\n";
    cout << "  Cannot open the input file " << input_name << "\n";
    return true;
  }
//
//  Read the header.
//
  error = pgmb_read_header ( input, xsize, ysize, maxg );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_READ: Fatal error!\n";
    cout << "  PGMB_READ_HEADER failed.\n";
    return true;
  }
//
//  Allocate storage for the data.
//
  *g = new unsigned char [ (*xsize) * (*ysize) ];
//
//  Read the data.
//
  error = pgmb_read_data ( input, *xsize, *ysize, *g );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_READ: Fatal error!\n";
    cout << "  PGMB_READ_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  input.close ( );

  return false;
}
//****************************************************************************80

bool pgmb_read_data ( ifstream &input, int xsize, int ysize, 
  unsigned char *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_READ_DATA reads the data in a binary portable gray map file.
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
//    portable gray map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, unsigned char *G, the array of XSIZE by YSIZE data values.
//
//    Output, bool PGMB_READ_DATA, is true if an error occurred.
//
{
  char c;
  bool error;
  int i;
  unsigned char *indexg;
  int j;

  indexg = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      input.read ( &c, 1 );
      *indexg = ( unsigned char ) c;
      indexg = indexg + 1;
      error = input.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PGMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading pixel (" 
          << i << ", " << j <<") \n";
        return true;
      }
    }
  }
  return false;
}
//****************************************************************************80

bool pgmb_read_header ( ifstream &input, int *xsize, int *ysize, 
  unsigned char *maxg )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_READ_HEADER reads the header of a binary portable gray map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the file containing the binary
//    portable gray map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, unsigned char *MAXG, the maximum gray value.
//
//    Output, bool PGMB_READ_HEADER, is true if an error occurred.
//
{
  int count;
  int fred;
  string line;
  int maxg2;
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
      cout << "PGMB_READ_HEADER - Fatal error!\n";
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

      if ( !s_eqi ( word, "P5" ) )
      {
        cout << "\n";
        cout << "PGMB_READ_HEADER - Fatal error.\n";
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
      line = rest;
      step = 3;
    }

    if ( step == 3 )
    {
      s_word_extract_first ( line, word, rest );

      if ( s_len_trim ( word ) <= 0 )
      {
        continue;
      }
      fred = atoi ( word.c_str ( ) );
      *maxg = ( unsigned char ) fred;
      line = rest;
      break;
    }
  }

  return false;
}
//****************************************************************************80

bool pgmb_to_pgma ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_TO_PGMA converts one PGMB file to PGMA format.
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
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input PGMB file.
//
//    Input, char *FILE_OUT_NAME, the name of the output PGMA file.
//
//    Output, bool HANDLE, is true if an error occurred.
//
{
  bool error;
  unsigned char *g;
  int *g2;
  unsigned char maxg;
  int xsize;
  int ysize;
//
//  Read the input file.
//
  error = pgmb_read ( file_in_name, &xsize, &ysize, &maxg, &g );

  if ( error ) 
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA: Fatal error!\n";
    cout << "  PGMB_READ failed.\n";
    return true;
  }
//
//  Check the data.
//
  error = pgmb_check_data ( xsize, ysize, maxg, g );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA: Fatal error!\n";
    cout << "  PGMB_CHECK_DATA reports bad data from the file.\n";

    delete [] g;
    return true;
  }
//
//  Convert the data.
//
  g2 = new int [ xsize *  ysize ];
  ucvec_to_i4vec ( xsize * ysize, g, g2 );
  delete [] g;
//
//  Write the output file.
//
  pgma_write ( file_out_name, xsize, ysize, g2 );

  delete [] g2;

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_TO_PGMA: Fatal error!\n";
    cout << "  PGMA_WRITE failed.\n";
    return true;
  }

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
//****************************************************************************80

void ucvec_to_i4vec ( int n, unsigned char *a, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    UCVEC_TO_I4VEC converts a vector of UNSIGNED CHAR's to an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items to convert.
//
//    Input, unsigned char *A, a pointer to a vector of unsigned chars.
//   
//    Input, int *B, a pointer to a vector of ints.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    *b = ( int ) *a;
    a++;
    b++;
  }

  return;
}
