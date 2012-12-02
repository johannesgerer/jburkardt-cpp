# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>

using namespace std;

# include "pgma_io.hpp"

//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
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

void pgma_example ( int xsize, int ysize, int *g )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_EXAMPLE sets up some data for an ASCII PGM file.
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
//    Output, int *G, the array of XSIZE by YSIZE gray values.
//
{
  int i;
  int *indexg;
  int j;
  int periods = 3;
  double PI = 3.14159265;
  double x;
  double y;

  indexg = g;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( ( double ) ( 2 * i ) ) / ( ( double ) ( ysize - 1 ) ) - 1.0;

    for ( j = 0; j < xsize; j++ )
    {
      x = ( 2.0 * PI * ( double ) ( periods * j ) ) / ( ( double ) ( xsize - 1 ) );

      *indexg = ( int ) ( 20.0 * ( sin ( x ) - y + 2 ) );
 
      indexg = indexg + 1;
    }
  }

  return;
}
//****************************************************************************80

void pgma_read ( string input_name, int &xsize, int &ysize, int &maxg, int **g )

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
//    22 July 2011
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_NAME, the name of the file.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &MAXG, the maximum gray value.
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
  numbytes = xsize * ysize * sizeof ( int );

  *g = new int[numbytes];
//
//  Read the data.
//
  pgma_read_data ( input, xsize, ysize, *g );
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

void pgma_read_header ( ifstream &input, int &xsize, int &ysize, int &maxg )

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
//    22 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the file.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &MAXG, the maximum gray value.
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
      xsize = atoi ( word.c_str ( ) );
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
      ysize = atoi ( word.c_str ( ) );
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
      maxg = atoi ( word.c_str ( ) );
      break;
    }

  }

  return;
}
//****************************************************************************80

void pgma_read_test ( string input_name )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_READ_TEST tests the ASCII PGM read routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_NAME, the name of the file.
//
{
  int *g;
  int maxg;
  int xsize;
  int ysize;

  g = NULL;
//
//  Read the data.
//
  pgma_read ( input_name, xsize, ysize, maxg, &g );

  cout << "\n";
  cout << "PGMA_READ_TEST:\n";
  cout << "  PGMA_READ was able to read \"" << input_name << "\".\n";
//
//  Check the data.
//
  pgma_check_data ( xsize, ysize, maxg, g );

  delete [] g;

  cout << "\n";
  cout << "PGMA_READ_TEST:\n";
  cout << "  PGMA_CHECK_DATA has approved the data from the file.\n";

  return;
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

void pgma_write_test ( string output_name )

//****************************************************************************80
//
//  Purpose:
//
//    PGMA_WRITE_TEST tests the ASCII PGM write routines.
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
{
  int *g;
  int maxg;
  int xsize;
  int ysize;

  xsize = 300;
  ysize = 300;
//
//  Allocate memory.
//
  g = new int[xsize * ysize];
//
//  Set the data.
//
  pgma_example ( xsize, ysize, g );
//
//  Write the data to the file.
//
  pgma_write ( output_name, xsize, ysize, g );

  delete [] g;

  return;
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
