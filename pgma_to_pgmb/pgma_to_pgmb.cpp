# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void i4vec_to_ucvec ( int n, int *a, unsigned char *b );
void pgma_check_data ( int xsize, int ysize, int maxg, int *g );
void pgma_read ( string input_name, int *xsize, int *ysize, int *maxg, int **g );
void pgma_read_data ( ifstream &input, int xsize, int ysize, int *g );
void pgma_read_header ( ifstream &input, int *xsize, int *ysize, int *maxg );
bool pgma_to_pgmb ( char *file_in_name, char *file_out_name );
bool pgmb_write ( string output_name, int xsize, int ysize, unsigned char *g );
bool pgmb_write_data ( ofstream &output, int xsize, int ysize, 
  unsigned char *g );
bool pgmb_write_header ( ofstream &output, int xsize, int ysize, 
  unsigned char maxg );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PGMA_TO_PGMB.
//
//  Discussion:
//
//    PGMA_TO_PGMB converts an ASCII portable gray map to a binary
//    portable gray map.
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
//    pgma_to_pgmb file.pgma file.pgmb
//
//  Parameters:
//
//    FILE.PGMA is the name of the input ASCII PGM file to be read.
//
//    FILE.PGMB is the name of the output binary PGM file to be created.
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
    cout << "PGMA_TO_PGMB:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Convert an ASCII PGM file to binary PGM format.\n";
  }
//
//  Get the specification for the input file.
//
  if ( argc < 2 ) 
  {
    cout << "\n";
    cout << "PGMA_TO_PGMB:\n";
    cout << "  Please enter the input PGMA file name:\n";
    
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
    cout << "PGMA_TO_PGMB:\n";
    cout << "  Please enter the output PGMB file name:\n";

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

  error = pgma_to_pgmb ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PGMA_TO_PGMB - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  if ( verbose )
  {
    cout << "\n";
    cout << "PGMA_TO_PGMB:\n";
    cout << "  Normal end of execution.\n";

    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void i4vec_to_ucvec ( int n, int *a, unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_TO_UCVEC converts an I4VEC into UNSIGNED CHAR's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items to convert.
//
//    Input, int *A, a pointer to a vector of ints.
//   
//    Input, unsigned char *B, a pointer to a vector of unsigned chars.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    *b = ( unsigned char ) *a;
    a++;
    b++;
  }
  return;
}
//****************************************************************************80

void pgma_check_data ( int xsize, int ysize, int maxg, int *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_CHECK_DATA checks the data for an ASCII PGM file.
//
//  Discussion:
//
//    XSIZE and YSIZE must be positive, the pointers must not be null,
//    and the data must be nonnegative and no greater than MAXG.
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
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXG, the maximum gray value.
//
//    Input, int *G, the array of XSIZE by YSIZE data values.
//
{
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cerr<< "\n";
    cerr << "PGMA_CHECK_DATA: Error!\n";
    cerr << "  XSIZE <= 0.\n";
    cerr << "  XSIZE = " << xsize << "\n";
    exit ( 1 );
  }

  if ( ysize <= 0 )
  {
    cerr << "\n";
    cerr << "PGMA_CHECK_DATA: Error!\n";
    cerr << "  YSIZE <= 0.\n";
    cerr << "  YSIZE = " << ysize << "\n";
    exit ( 1 );
  }

  if ( g == NULL )
  {
    cerr << "\n";
    cerr << "PGMA_CHECK_DATA: Error!\n";
    cerr << "  Null pointer to g.\n";
    exit ( 1 );
  }

  index = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *index < 0 )
      {
        cerr << "\n";
        cerr << "PGMA_CHECK_DATA - Fatal error!\n";
        cerr << "  Negative data.\n";
        cerr << "  G(" << i << "," << j << ")=" << *index << "\n";
        exit ( 1 );
      }
      else if ( maxg < *index )
      {
        cerr << "\n";
        cerr << "PGMA_CHECK_DATA - Fatal error!\n";
        cerr << "  Data exceeds MAXG = " << maxg << "\n";
        cerr << "  G(" << i << "," << j << ")=" << *index << "\n";
        exit ( 1 );
      }
      index = index + 1;
    }
  } 
  return;
}
//****************************************************************************80

void pgma_read ( string input_name, int *xsize, int *ysize, int *maxg,
  int **g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_READ reads the header and data from an ASCII PGM file.
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
//    Input, string INPUT_NAME, the name of the file.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXG, the maximum gray value.
//
//    Output, int **G, the array of XSIZE by YSIZE data values.
//
{
  ifstream input;
  int numbytes;

  input.open ( input_name.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "PGMA_READ - Fatal error!\n";
    cerr << "  Cannot open the input file \"" << input_name << "\".\n";
    exit ( 1 );
  }
//
//  Read the header.
//
  pgma_read_header ( input, xsize, ysize, maxg );
//
//  Allocate storage for the data.
//
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *g = new int[numbytes];
//
//  Read the data.
//
  pgma_read_data ( input, *xsize, *ysize, *g );
//
//  Close the file.
//
  input.close ( );

  return;
}
//****************************************************************************80

void pgma_read_data ( ifstream &input, int xsize, int ysize, int *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_READ_DATA reads the data in an ASCII PGM file.
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
//    Input, ifstream &INPUT, a pointer to the file containing the data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *G, the array of XSIZE by YSIZE data values.
//
{
  int i;
  int j;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      input >> *g;
      if ( input.eof ( ) )
      {
        exit ( 1 );
      }
      g = g + 1;
    }
  }

  return;
}
//****************************************************************************80

void pgma_read_header ( ifstream &input, int *xsize, int *ysize, int *maxg )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_READ_HEADER reads the header of an ASCII PGM file.
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
//    Input, ifstream &INPUT, a pointer to the file.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXG, the maximum gray value.
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
      cerr << "\n";
      cerr << "PGMA_READ_HEADER - Fatal error!\n";
      cerr << "  End of file.\n";
      exit ( 1 );
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    if ( step == 0 )
    {
      s_word_extract_first ( line, word, rest );

      if ( s_len_trim ( word ) == 0 )
      {
        continue;
      }
      line = rest;

      if ( ( word[0] != 'P' && word[0] != 'p' ) || 
             word[1] != '2' )
      {
        cerr << "\n";
        cerr << "PGMA_READ_HEADER - Fatal error.\n";
        cerr << "  Bad magic number = \"" << word << "\".\n";
        exit ( 1 );
      }
      step = 1;
    }

    if ( step == 1 )
    {
      s_word_extract_first ( line, word, rest );

      if ( s_len_trim ( word ) == 0 )
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

      if ( s_len_trim ( word ) == 0 )
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

      if ( s_len_trim ( word ) == 0 )
      {
        continue;
      }
      *maxg = atoi ( word.c_str ( ) );
      break;
    }

  }

  return;
}
//****************************************************************************80

bool pgma_to_pgmb ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_TO_PGMB converts an ASCII PGM file to binary PGM format.
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
//    Input, char *FILE_IN_NAME, the name of the ASCII PGM file to be read.
//
//    Input, char *FILE_OUT_NAME, the name of the binary PGM file to be created.
//   
//    Output, bool PGMA_TO_PGMB, is true if an error occurred.
//
{
  int *g;
  unsigned char *g2;
  bool error;
  int maxg;
  int xsize;
  int ysize;
//
//  Read the input file.
//
  pgma_read ( file_in_name, &xsize, &ysize, &maxg, &g );

  if ( error )
  {
    cout << "\n";
    cout << "PGMA_TO_PGMB: Fatal error!\n";
    cout << "  PGMA_READ failed.\n";
    return true;
  }
//
//  Check the data.
//
  pgma_check_data ( xsize, ysize, maxg, g );
//
//  Copy the data into unsigned char's.
//
  g2 = new unsigned char [ xsize * ysize ];
  i4vec_to_ucvec ( xsize * ysize, g, g2 );
  delete [] g;
//
//  Write the output file.
//
  error = pgmb_write ( file_out_name, xsize, ysize, g2 );

  delete [] g2;

  if ( error )
  {
    cout << "\n";
    cout << "PGMA_TO_PGMB: Fatal error!\n";
    cout << "  PGMB_WRITE failed.\n";
    return true;
  }

  return false;
}
//****************************************************************************80

bool pgmb_write ( string output_name, int xsize, int ysize, unsigned char *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_WRITE writes the header and data for a binary portable gray map file.
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
//    Input, string OUTPUT_NAME, the name of the file to contain the binary
//    portable gray map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *G, the array of XSIZE by YSIZE data values.
//
//    Output, bool PGMB_WRITE, is true if an error occurred.
//
{
  bool error;
  ofstream output;
  int i;
  unsigned char *indexg;
  int j;
  unsigned char maxg;
//
//  Determine the maximum gray value.
//
  maxg = 0;
  indexg = g;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      if ( maxg < *indexg )
      {
        maxg = *indexg;
      }
      indexg = indexg + 1;
    }
  }
//
//  Open the file.
//
  output.open ( output_name.c_str ( ), ios::binary );

  if ( !output )
  {
    cout << "\n";
    cout << "PGMB_WRITE: Fatal error!\n";
    cout << "  Cannot open the output file " << output_name << "\n";
    return true;
  }
//
//  Write the header.
//
  error = pgmb_write_header ( output, xsize, ysize, maxg );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_WRITE: Fatal error!\n";
    cout << "  PGMB_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = pgmb_write_data ( output, xsize, ysize, g );

  if ( error )
  {
    cout << "\n";
    cout << "PGMB_WRITE: Fatal error!\n";
    cout << "  PGMB_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  output.close ( );

  return false;
}
//****************************************************************************80

bool pgmb_write_data ( ofstream &output, int xsize, int ysize, 
  unsigned char *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_WRITE_DATA writes the data for a binary portable gray map file.
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
//    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
//    portable gray map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char *G, the array of XSIZE by YSIZE data values.
//
//    Output, bool PGMB_WRITE_DATA, is true if an error occurred.
//
{
  int i;
  unsigned char *indexg;
  int j;

  indexg = g;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      output << *indexg;
      indexg = indexg + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool pgmb_write_header ( ofstream &output, int xsize, int ysize, 
  unsigned char maxg )

//****************************************************************************80
//
//  Purpose:
//
//    PGMB_WRITE_HEADER writes the header of a binary portable gray map file.
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
//    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
//    portable gray map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char MAXG, the maximum gray value.
//
//    Output, bool PGMB_WRITE_HEADER, is true if an error occurred.
//
{
  output << "P5" << " "
           << xsize << " " 
           << ysize << " " 
           << ( int ) maxg << "\n";

  return false;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM2 returns the length of a string to the last nonblank.
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
//    Output, int S_LEN_TRIM2, the length of the string to the last nonblank.
//    If S_LEN_TRIM2 is 0, then the string is entirely blank.
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
//    S_WORD_EXTRACT_FIRST2 extracts the first word from a string.
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
