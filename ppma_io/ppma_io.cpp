# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "ppma_io.hpp"

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
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
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

bool ppma_check_data ( int xsize, int ysize, int rgb_max, int *r,
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_CHECK_DATA checks the data for an ASCII portable pixel map file.
//
//  Discussion:
//
//    XSIZE and YSIZE must be positive, the pointers must not be null,
//    and the data must be nonnegative and no greater than RGB_MAX.
//
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int RGB_MAX, the maximum RGB value.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_CHECK_DATA, is
//    true, if an error was detected, or
//    false, if the data was legal.
//
{
  char c;
  int i;
  int *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Fatal error!\n";
    cout << "  xsize <= 0.\n";
    cout << "  xsize = " << xsize << "\n";
    return true;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Fatal error!\n";
    cout << "  ysize <= 0.\n";
    cout << "  ysize = " << ysize << "\n";
    return true;
  }

  if ( r == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Fatal error!\n";
    cout << "  Null pointer to R.\n";
    return true;
  }

  if ( g == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Fatal error!\n";
    cout << "  Null pointer to G.\n";
    return true;
  }

  if ( b == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Fatal error!\n";
    cout << "  Null pointer to B.\n";
    return true;
  }

  for ( k = 0; k < 3; k++ )
  {

    if ( k == 0 )
    {
      index = r;
      c = 'R';
    }
    else if ( k == 1 )
    {
      index = g;
      c = 'G';
    }
    else if ( k == 2 )
    {
      index = b;
      c = 'B';
    }

    for ( j = 0; j < ysize; j++ )
    {
      for ( i = 0; i < xsize; i++ )
      {
        if ( *index < 0 )
        {
          cout << "\n";
          cout << "PPMA_CHECK_DATA - Fatal error!\n";
          cout << "  Negative data.\n";
          cout << "  " << c << "(" << i << "," << j << ")=" << *index << "\n";
          return true;
        }
        else if ( rgb_max < *index )
        {
          cout << "\n";
          cout << "PPMA_CHECK_DATA - Fatal error!\n";
          cout << "  Data exceeds RGB_MAX = " << rgb_max << "\n";
          cout << "  " << c << "(" << i << "," << j << ")=" << *index << "\n";
          return true;
        }

        index = index + 1;
      }
    } 
  }

  return false;
}
//****************************************************************************80

bool ppma_example ( int xsize, int ysize, int *r, int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_EXAMPLE sets up some RGB data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE RGB values.
//
//    Output, bool PPMA_EXAMPLE, is
//    false, if no error occurred,
//    true, if an error occurred.
//
{
  int *b_index;
  float f1;
  float f2;
  float f3;
  int *g_index;
  int i;
  int j;
  int *r_index;
  float x;
  float y;

  r_index = r;
  g_index = g;
  b_index = b;

  for ( i = 0; i < ysize; i++ )
  {
    y = ( float ) ( ysize + 1 - i ) / ( float ) ( ysize - 1 );
    for ( j = 0; j < xsize; j++ )
    {
      x = ( float ) ( j ) / ( float ) ( xsize - 1 );

      f1 = 4.0 * ( x - 0.5 ) * ( x - 0.5 );
      f2 = sin ( 3.14159265 * x );
      f3 = x;

      if ( y <= f1 )
      {
        *r_index = ( int ) ( 255.0 * f1 );
      }
      else
      {
        *r_index = 50;
      }

      if ( y <= f2 )
      {
        *g_index = ( int ) ( 255.0 * f2 );
      }
      else
      {
        *g_index = 150;
      }

      if ( y <= f3 )
      {
        *b_index = ( int ) ( 255.0 * f3 );
      }
      else
      {
        *b_index = 250;
      }

      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool ppma_read ( string input_name, int &xsize, int &ysize, int &rgb_max,
  int **r, int **g, int **b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_READ reads the header and data from an ASCII portable pixel map file.
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
//    Input, string INPUT_NAME, the name of the file containing the ASCII
//    portable pixel map data.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &RGB_MAX, the maximum RGB value.
//
//    Output, int **R, **G, **B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_READ, is
//    true, if an error was detected, or
//    false, if the file was read.
//
{
  bool error;
  ifstream input;
  int numbytes;

  input.open ( input_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_name << "\".\n";
    return true;
  }
//
//  Read the header.
//
  error = ppma_read_header ( input, xsize, ysize, rgb_max );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  PPMA_READ_HEADER failed.\n";
    return true;
  }
//
//  Allocate storage for the data.
//
  numbytes = xsize * ysize * sizeof ( int );

  *r = new int[numbytes];
  *g = new int[numbytes];
  *b = new int[numbytes];
//
//  Read the data.
//
  error = ppma_read_data ( input, xsize, ysize, *r, *g, *b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  PPMA_READ_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  input.close ( );

  return false;
}
//****************************************************************************80

bool ppma_read_data ( ifstream &input, int xsize, int ysize, int *r,
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_READ_DATA reads the data in an ASCII portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &INPUT, a pointer to the file containing the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Output, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_READ_DATA, is
//    true, if an error was detected, or
//    false, if the data was read.
//
{
  int i;
  int j;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      input >> *r;
      if ( input.eof() )
      {
        return true;
      }
      r = r + 1;

      input >> *g;
      if ( input.eof() )
      {
        return true;
      }
      g = g + 1;

      input >> *b;
      if ( input.eof() )
      {
        return true;
      }
      b = b + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool ppma_read_header ( ifstream &input, int &xsize, int &ysize, 
  int &rgb_max )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_READ_HEADER reads the header of an ASCII portable pixel map file.
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
//    Input, ifstream &INPUT, a pointer to the file containing the ASCII
//    portable pixel map data.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &RGB_MAX, the maximum RGB value.
//
//    Output, bool PPMA_READ_HEADER, is
//    true, if an error was detected, or
//    false, if the header was read.
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
      cout << "PPMA_READ_HEADER - Fatal error!\n";
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

      if ( !s_eqi ( word, "P3" ) )
      {
        cout << "\n";
        cout << "PPMA_READ_HEADER - Fatal error.\n";
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
      xsize = atoi ( word.c_str ( ) );
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
      ysize = atoi ( word.c_str ( ) );
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
      rgb_max = atoi ( word.c_str ( ) );
      line = rest;
      break;
    }

  }

  return false;
}
//****************************************************************************80

bool ppma_read_test ( string input_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_READ_TEST tests the ASCII portable pixel map read routines.
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
//    Input, string INPUT_NAME, the name of the file containing the ASCII
//    portable pixel map data.
//
//    Output, bool PPMA_READ_TEST, is
//    true, if an error was detected, or
//    false, if the test was carried out.
//
{
  int *b;
  bool error;
  int *g;
  int *r;
  int rgb_max;
  int xsize;
  int ysize;

  r = NULL;
  g = NULL;
  b = NULL;
//
//  Read the data.
//
  error = ppma_read ( input_name, xsize, ysize, rgb_max, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ_TEST - Fatal error!\n";
    cout << "  PPMA_READ failed.\n";

    delete [] r;
    delete [] g;
    delete [] b;

    return true;
  }

  cout << "\n";
  cout << "PPMA_READ_TEST:\n";
  cout << "  PPMA_READ was able to read \"" << input_name << "\".\n";
//
//  Check the data.
//
  error = ppma_check_data ( xsize, ysize, rgb_max, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_READ_TEST - Fatal error!\n";
    cout << "  PPMA_CHECK_DATA reports bad data in the file.\n";
    return true;
  }

  cout << "\n";
  cout << "PPMA_READ_TEST:\n";
  cout << "  PPMA_CHECK_DATA has approved the data from the file.\n";

  return false;
}
//****************************************************************************80

bool ppma_write ( string output_name, int xsize, int ysize, int *r, 
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.
// 
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    28 February 2003
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE, is
//    true, if an error was detected, or
//    false, if the file was written.
//
{
  int *b_index;
  bool error;
  ofstream output;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
//
//  Open the output file.
//
  output.open ( output_name.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_name << "\".\n";
    return true;
  }
//
//  Compute the maximum.
//
  rgb_max = 0;
  r_index = r;
  g_index = g;
  b_index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( rgb_max < *r_index )
      {
        rgb_max = *r_index;
      }
      r_index = r_index + 1;

      if ( rgb_max < *g_index )
      {
        rgb_max = *g_index;
      }
      g_index = g_index + 1;

      if ( rgb_max < *b_index )
      {
        rgb_max = *b_index;
      }
      b_index = b_index + 1;
    }
  }
//
//  Write the header.
//
  error = ppma_write_header ( output, output_name, xsize, ysize, rgb_max );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = ppma_write_data ( output, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  output.close ( );

  return false;
}
//****************************************************************************80

bool ppma_write_data ( ofstream &output, int xsize, int ysize, int *r,
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE_DATA, is
//    true, if an error was detected, or
//    false, if the data was written.
//
{
  int *b_index;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_num;

  r_index = r;
  g_index = g;
  b_index = b;
  rgb_num = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      output << *r_index << " " << *g_index << " " << *b_index;
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        output << "\n";
      }
      else
      {
        output << " ";
      }
    }
  }
  return false;
}
//****************************************************************************80

bool ppma_write_header ( ofstream &output, string output_name, int xsize, 
  int ysize, int rgb_max )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, string OUTPUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int RGB_MAX, the maximum RGB value.
//
//    Output, bool PPMA_WRITE_HEADER, is
//    true, if an error was detected, or
//    false, if the header was written.
//
{
  output << "P3\n";
  output << "# " << output_name << " created by PPMA_IO::PPMA_WRITE.C.\n";
  output << xsize << "  " << ysize << "\n";
  output << rgb_max << "\n";

  return false;
}
//****************************************************************************80

bool ppma_write_test ( string output_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_TEST tests the ASCII portable pixel map write routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Output, bool PPMA_WRITE_TEST, equals
//    true, if the test could not be carried out,
//    false, if the test was carried out.
//
{
  int *b;
  bool error;
  int *g;
  int *r;
  int xsize;
  int ysize;

  xsize = 300;
  ysize = 300;
//
//  Allocate memory.
//
  r = new int[xsize * ysize];
  g = new int[xsize * ysize];
  b = new int[xsize * ysize];
//
//  Set the data.
//
  error = ppma_example ( xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE_TEST - Fatal error!\n";
    cout << "  PPMA_EXAMPLE failed.\n";
    return true;
  }
//
//  Write the data to the file.
//
  error = ppma_write ( output_name, xsize, ysize, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;
 
  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE_TEST - Fatal error!\n";
    cout << "  PPMA_WRITE failed.\n";
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
