# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void pbma_check_data ( int xsize, int ysize, int *b );
void pbma_read ( string file_in_name, int *xsize, int *ysize, int **b );
void pbma_read_data ( ifstream &file_in, int xsize, int ysize, int *b );
void pbma_read_header ( ifstream &file_in, int *xsize, int *ysize );
void pbma_to_pbmb ( string input_file, string output_file );
bool pbmb_write ( string output_name, int xsize, int ysize, int *barray );
bool pbmb_write_data ( ofstream &output, int xsize, int ysize, int *barray );
bool pbmb_write_header ( ofstream &output, int xsize, int ysize );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PBMA_TO_PBMB.
//
//  Discussion:
//
//    PBMA_TO_PBMB converts an ASCII portable bit map to a binary
//    portable bit map.
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
//  Usage:
//
//    pbma_to_pbmb file.pbma file.pbmb
//
//  Parameters:
//
//    FILE.PBMA is the name of the input ASCII PBM file to be read.
//
//    FILE.PBMB is the name of the output binary PBM file to be created.
//
{
  bool error;
  string input_file;
  string output_file;
  bool verbose = false;

  if ( verbose )
  {
    timestamp ( );
    cout << "\n";
    cout << "PBMA_TO_PBMB:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Convert an ASCII PBM file to binary PBM format.\n";
  }
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PBMA_TO_PBMB:\n";
    cout << "  Please enter the input ASCII PBM file name:\n";

    cin >> input_file;

    error = cin.eof ( );

    if ( error )
    {
      exit ( 1 );
    }
  }
  else
  {
    input_file = argv[1];
  }
//
//  Get the specification for the output file.
//
  if ( argc < 3 )
  {
    cout << "\n";
    cout << "PBMA_TO_PBMB:\n";
    cout << "  Please enter the output binary PBM file name:\n";

    cin >> output_file;

    error = cin.eof ( );

    if ( error )
    {
      exit ( 1 );
    }
  }
  else
  {
    output_file = argv[2];
  }

  pbma_to_pbmb ( input_file, output_file );

  if ( verbose )
  {
    cout << "\n";
    cout << "PBMA_TO_PBMB:\n";
    cout << "  Normal end of execution.\n";

    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void pbma_check_data ( int xsize, int ysize, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_CHECK_DATA checks the data for an ASCII PBM file.
//
//  Discussion:
//
//    XSIZE and YSIZE must be positive, the pointers must not be null,
//    and the data must be 0 or 1.
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
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *B, the array of XSIZE by YSIZE data values.
//
{
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cerr << "\n";
    cerr << "PBMA_CHECK_DATA: Error!\n";
    cerr << "  XSIZE <= 0.\n";
    cerr << "  XSIZE = " << xsize << "\n";
    exit ( 1 );
  }

  if ( ysize <= 0 )
  {
    cerr << "\n";
    cerr << "PBMA_CHECK_DATA: Error!\n";
    cerr << "  YSIZE <= 0.\n";
    cerr << "  YSIZE = " << ysize << "\n";
    exit ( 1 );
  }

  if ( b == NULL )
  {
    cerr << "\n";
    cerr << "PBMA_CHECK_DATA: Error!\n";
    cerr << "  Null pointer to B.\n";
    exit ( 1 );
  }

  index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( *index < 0 )
      {
        cerr << "\n";
        cerr << "PBMA_CHECK_DATA - Fatal error!\n";
        cerr << "  Negative data.\n";
        cerr << "  B(" << i << "," << j << ")=" << *index << "\n";
        exit ( 1 );
      }
      else if ( 1 < *index )
      {
        cerr << "\n";
        cerr << "PBMA_CHECK_DATA - Fatal error!\n";
        cerr << "  Data exceeds 1\n";
        cerr << "  B(" << i << "," << j << ")=" << *index << "\n";
        exit ( 1 );
      }

      index = index + 1;
    }
  }
  return;
}
//****************************************************************************80

void pbma_read ( string file_in_name, int *xsize, int *ysize, int **b )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_READ reads the header and data from an ASCII PBM file.
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
//    Input, string FILE_IN_NAME, the name of the file.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int **B, the array of XSIZE by YSIZE data values.
//
{
  ifstream file_in;
  int numbytes;

  file_in.open ( file_in_name.c_str ( ) );

  if ( !file_in )
  {
    cerr << "\n";
    cerr << "PBMA_READ - Fatal error!\n";
    cerr << "  Cannot open the input file \"" << file_in_name << "\".\n";
    exit ( 1 );
  }
//
//  Read the header.
//
  pbma_read_header ( file_in, xsize, ysize );
//
//  Allocate storage for the data.
//
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *b = new int[numbytes];
//
//  Read the data.
//
  pbma_read_data ( file_in, *xsize, *ysize, *b );
//
//  Close the file.
//
  file_in.close ( );

  return;
}
//****************************************************************************80

void pbma_read_data ( ifstream &file_in, int xsize, int ysize, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_READ_DATA reads the data in an ASCII PBM file.
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
//    Input, ifstream &FILE_IN, a pointer to the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *B, the array of XSIZE by YSIZE data values.
//
{
  int i;
  int j;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      file_in >> *b;
      if ( file_in.eof ( ) )
      {
        exit ( 1 );
      }
      b = b + 1;
    }
  }
  return;
}
//****************************************************************************80

void pbma_read_header ( ifstream &file_in, int *xsize, int *ysize )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_READ_HEADER reads the header of an ASCII PBM file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_IN, a pointer to the file.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
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
    getline ( file_in, line );

    if ( file_in.eof ( ) )
    {
      cerr << "\n";
      cerr << "PBMA_READ_HEADER - Fatal error!\n";
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
             word[1] != '1' )
      {
        cerr << "\n";
        cerr << "PBMA_READ_HEADER - Fatal error.\n";
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
      break;
    }
  }
  return;
}
//****************************************************************************80

void pbma_to_pbmb ( string input_file, string output_file )

//****************************************************************************80
//
//  Purpose:
//
//    PBMA_TO_PBMB converts an ASCII PBM file to binary PBM format.
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
//    Input, string INPUT_FILE, the name of the ASCII PBM file to be read.
//
//    Input, string OUTPUT_FILE, the name of the binary PBM file to be created.
//
{
  int *b;
  bool error;
  int xsize;
  int ysize;
//
//  Read the input file.
//
  pbma_read ( input_file, &xsize, &ysize, &b );
//
//  Check the data.
//
  pbma_check_data ( xsize, ysize, b );
//
//  Write the output file.
//
  error = pbmb_write ( output_file, xsize, ysize, b );

  delete [] b;

  if ( error )
  {
    cout << "\n";
    cout << "PBMA_TO_PBMB: Fatal error!\n";
    cout << "  PBMB_WRITE failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

bool pbmb_write ( string output_name, int xsize, int ysize, int *barray )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_WRITE writes the header and data for a binary portable bit map file.
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
//    portable bit map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *BARRAY, the array of XSIZE by YSIZE data values.
//
//    Output, bool PBMB_WRITE, is true if an error occurred.
//
{
  bool error;
  ofstream output;

  output.open ( output_name.c_str ( ), ios::binary );

  if ( !output )
  {
    cout << "\n";
    cout << "PBMB_WRITE: Fatal error!\n";
    cout << "  Cannot open the output file " << output_name << "\n";
    return true;
  }
//
//  Write the header.
//
  error = pbmb_write_header ( output, xsize, ysize );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_WRITE: Fatal error!\n";
    cout << "  PBMB_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = pbmb_write_data ( output, xsize, ysize, barray );

  if ( error )
  {
    cout << "\n";
    cout << "PBMB_WRITE: Fatal error!\n";
    cout << "  PBMB_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  output.close ( );

  return false;
}
//****************************************************************************80

bool pbmb_write_data ( ofstream &output, int xsize, int ysize, int *barray )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_WRITE_DATA writes the data for a binary portable bit map file.
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
//    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
//    portable bit map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *BARRAY, the array of XSIZE by YSIZE data values.
//
//    Output, bool PBMB_WRITE_DATA, is true if an error occurred.
//
{
  int bit;
  unsigned char c;
  int i;
  int *indexb;
  int j;
  int k;

  indexb = barray;
  c = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      k = 7 - i%8;
      bit = (*indexb)%2;
      c = c | ( bit << k );

      indexb = indexb + 1;

      if ( (i+1)%8 == 0 || i == ( xsize - 1 ) )
      {
        output << c;
        c = 0;
      }
    }
  }
  return false;
}
//****************************************************************************80

bool pbmb_write_header ( ofstream &output, int xsize, int ysize )

//****************************************************************************80
//
//  Purpose:
//
//    PBMB_WRITE_HEADER writes the header of a binary portable bit map file.
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
//    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
//    portable bit map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, bool PBMB_WRITE_HEADER, is true if an error occurred.
//
{
  output << "P4" << " "
           << xsize << " "
           << ysize << "\n";
 
  return false;
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
