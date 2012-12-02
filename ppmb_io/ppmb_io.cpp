# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "ppmb_io.hpp"

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

bool ppmb_check_data ( int xsize, int ysize, int maxrgb, unsigned char *r,
  unsigned char *g, unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_CHECK_DATA checks the data for a binary portable pixel map file.
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
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXRGB, the maximum RGB value.
//
//    Input, unsigned char *R, *G, *B, the arrays of XSIZE by YSIZE 
//    data values.
//
//    Output, bool PPM_CHECK_DATA, is true if an error occurred.
//
{
  int i;
  unsigned char *index;
  int j;
  int k;

  if ( xsize <= 0 )
  {
    cout << "\n";
    cout << "PPMB_CHECK_DATA - Error!\n";
    cout << "  XSIZE <= 0\n";
    cout << "  XSIZE = " << xsize << "\n";
    return 1;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PPMB_CHECK_DATA - Error!\n";
    cout << "  YSIZE <= 0\n";
    cout << "  YSIZE = " << ysize << "\n";
    return 1;
  }

  if ( r == NULL || g == NULL || b == NULL )
  {
    cout << "\n";
    cout << "PPMB_CHECK_DATA - Error!\n";
    cout << "  Null pointer to data.\n";
    return 1;
  }

  for ( k = 0; k < 3; k++ )
  {
    if ( k == 0 )
    {
      index = r;
    }
    else if ( k == 1 )
    {
      index = g;
    }
    else if ( k == 2 )
    {
      index = b;
    }

    for ( j = 0; j < ysize; j++ )
    {
      for ( i = 0; i < xsize; i++ )
      {
        if ( maxrgb < *index )
        {
          cout << "\n";
          cout << "PPMB_CHECK_DATA - Error!\n";
          cout << "  Pixel color value exceeds MAXRGB.\n";
          cout << "  MAXRGB = " << maxrgb << "\n";
          if ( k == 0 )
          {
            cout << "  R(" << i << "," << j << ") = " << *index << "\n";
          }
          else if ( k == 1 )
          {
            cout << "  G(" << i << "," << j << ") = " << *index << "\n";
          }
          else if ( k == 2 )
          {
            cout << "  B(" << i << "," << j << ") = " << *index << "\n";
          }
          return true;
        }
        index = index + 1;
      }
    } 
  }

  return 0;
}
//****************************************************************************80

bool ppmb_example ( int xsize, int ysize, unsigned char *r, unsigned char *g, 
  unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_EXAMPLE sets up some RGB data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int xsize, ysize, the number of rows and columns of data.
//
//    Output, unsigned char *r, *g, *b, the arrays of xsize by ysize RGB values.
//
//    Output, bool PPMA_EXAMPLE, is true if an error occurred.
//
{
  float f1;
  float f2;
  float f3;
  int i;
  unsigned char *indexr;
  unsigned char *indexg;
  unsigned char *indexb;
  int j;
  float x;
  float y;

  indexr = r;
  indexg = g;
  indexb = b;

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
        *indexr = ( unsigned char ) ( 255.0 * f1 );
      }
      else
      {
        *indexr = 50;
      }

      if ( y <= f2 )
      {
        *indexg = ( unsigned char ) ( 255.0 * f2 );
      }
      else
      {
        *indexg = 150;
      }

      if ( y <= f3 )
      {
        *indexb = ( unsigned char ) ( 255.0 * f3 );
      }
      else
      {
        *indexb = 250;
      }

      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool ppmb_read ( string input_name, int &xsize, int &ysize, int &maxrgb,
  unsigned char **r, unsigned char **g, unsigned char **b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_READ reads the header and data from a binary portable pixel map file.
//
//  Discussion:
//
//    Thanks to Jonas Schwertfeger for pointing out that, especially on 
//    Microsoft Windows systems, a binary file needs to be opened as a binary 
//    file!
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
//    Input, string INPUT_NAME, the name of the file containing the binary
//    portable pixel map data.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &MAXRGB, the maximum RGB value.
//
//    Output, unsigned char **R, **G, **B, the arrays of XSIZE by YSIZE 
//    data values.
//
//    Output, bool PPMB_READ, is true if there was an error.
//
{
  bool error;
  ifstream input;

  input.open ( input_name.c_str ( ), ios::binary );

  if ( !input )
  {
    cout << "\n";
    cout << "PPMB_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << input_name << "\".\n";
    return true;
  }
//
//  Read the header.
//
  error = ppmb_read_header ( input, xsize, ysize, maxrgb );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_READ - Fatal error!\n";
    cout << "  PPMB_READ_HEADER failed.\n";
    return true;
  }
//
//  Allocate storage for the data.
//
  *r = new unsigned char[xsize * ysize];
  *g = new unsigned char[xsize * ysize];
  *b = new unsigned char[xsize * ysize];
//
//  Read the data.
//
  error = ppmb_read_data ( input, xsize, ysize, *r, *g, *b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_READ - Fatal error!\n";
    cout << "  PPMB_READ_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  input.close ( );

  return false;
}
//****************************************************************************80

bool ppmb_read_data ( ifstream &input, int xsize, int ysize, 
  unsigned char *r, unsigned char *g, unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_READ_DATA reads the data in a binary portable pixel map file.
//
//  Discussion:
//
//    If the ordinary ">>" operator is used to input the data, then data that
//    happens to look like new lines or other white space is skipped.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &input, a pointer to the file containing the binary
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char *R, *G, *B, the arrays of XSIZE by YSIZE 
//    data values.
//
//    Output, bool PPMB_READ_DATA, is true if an error occurred.
//
{
  char c;
  bool error;
  int i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;
  int k;

  indexr = r;
  indexg = g;
  indexb = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      input.read ( &c, 1 );
      *indexr = ( unsigned char ) c;
      indexr = indexr + 1;
      error = input.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PPMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading R byte of pixel (" 
          << i << ", " << j <<") \n";
        return true;
      }

      input.read ( &c, 1 );
      *indexg = ( unsigned char ) c;
      indexg = indexg + 1;
      error = input.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PPMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading G byte of pixel (" 
          << i << ", " << j <<") \n";
        return true;
      }

      input.read ( &c, 1 );
      *indexb = ( unsigned char ) c;
      indexb = indexb + 1;
      error = input.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PPMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading B byte of pixel (" 
          << i << ", " << j <<") \n";
        return true;
      }

    }
  }
  return false;
}
//****************************************************************************80

bool ppmb_read_header ( ifstream &input, int &xsize, int &ysize, int &maxrgb )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_READ_HEADER reads the header of a binary portable pixel map file.
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
//    Input, ifstream &INPUT, a pointer to the file containing the binary
//    portable pixel map data.
//
//    Output, int &XSIZE, &YSIZE, the number of rows and columns of data.
//
//    Output, int &MAXRGB, the maximum RGB value.
//
//    Output, bool PPMB_READ_HEADER, is true if an error occurred.
//
{
  int count;
  int fred;
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
      cout << "PPMB_READ_HEADER - Fatal error!\n";
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

      if ( !s_eqi ( word, "P6" ) )
      {
        cout << "\n";
        cout << "PPMB_READ_HEADER - Fatal error.\n";
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
      fred = atoi ( word.c_str ( ) );
      maxrgb = ( unsigned char ) fred;
      line = rest;
      break;
    }

  }

  return false;
}
//****************************************************************************80

bool ppmb_read_test ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_READ_TEST tests the binary portable pixel map read routines.
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
//    Input, string FILE_NAME, the name of the file containing the binary
//    portable pixel map data.
//
//    Output, bool PPMB_READ_TEST, is true if an error occurred.
//
{
  unsigned char *b;
  bool error;
  unsigned char *g;
  int maxrgb;
  unsigned char *r;
  int xsize;
  int ysize;
//
//  Read the data.
//
  error = ppmb_read ( file_name, xsize, ysize, maxrgb, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_READ_TEST: Fatal error!\n";
    cout << "  PPMB_READ failed.\n";
    delete [] r;
    delete [] g;
    delete [] b;
    return true;
  }
//
//  Check the data.
//
  error = ppmb_check_data ( xsize, ysize, maxrgb, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;

  if ( error )
  {
    cout << "\n";
    cout << "  PPMB_READ_TEST - Warning\n";
    cout << "    PPM_CHECK data reports bad data from the file.\n";
    return true;
  }

  cout << "\n";
  cout << "  PPMB_READ_TEST:\n";
  cout << "  PPM_CHECK_DATA passes the data from the file.\n";

  return false;
}
//****************************************************************************80

bool ppmb_write ( string output_name, int xsize, int ysize, 
  unsigned char *r, unsigned char *g, unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_WRITE writes the header and data for a binary portable pixel map file.
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
//    11 April 2003
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_NAME, the name of the file to contain the binary
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char *R, *G, *B, the arrays of XSIZE by YSIZE 
//    data values.
//
//    Output, bool PPMB_WRITE, is true if an error occurred.
//
{
  bool error;
  ofstream output;
  int i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;
  int maxrgb;
//
//  Open the output file.
//
  output.open ( output_name.c_str ( ), ios::binary );

  if ( !output )
  {
    cout << "\n";
    cout << "PPMB_WRITE: Fatal error!\n";
    cout << "  Cannot open the output file " << output_name << ".\n";
    return true;
  }
//
//  Compute the maximum.
//
  maxrgb = 0;
  indexr = r;
  indexg = g;
  indexb = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      maxrgb = i4_max ( maxrgb, *indexr );
      maxrgb = i4_max ( maxrgb, *indexg );
      maxrgb = i4_max ( maxrgb, *indexb );
      indexr = indexr + 1;
      indexg = indexg + 1;
      indexb = indexb + 1;
    }
  }
//
//  Write the header.
//
  error = ppmb_write_header ( output, xsize, ysize, maxrgb );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_WRITE: Fatal error!\n";
    cout << "  PPMB_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = ppmb_write_data ( output, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_WRITE: Fatal error!\n";
    cout << "  PPMB_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  output.close ( );

  return false;
}
//****************************************************************************80

bool ppmb_write_data ( ofstream &output, int xsize, int ysize, 
  unsigned char *r, unsigned char *g, unsigned char *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_WRITE_DATA writes the data for a binary portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream OUTPUT, a pointer to the file to contain the binary
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, unsigned char *R, *G, *B, the arrays of XSIZE by YSIZE 
//    data values.
//
//    Output, bool PPMB_WRITE_DATA, is true if an error occurred.
//
{
  int  i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int  j;

  indexr = r;
  indexg = g;
  indexb = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      output << *indexr;
      indexr = indexr + 1;
      output << *indexg;
      indexg = indexg + 1;
      output << *indexb;
      indexb = indexb + 1;
    }
  }
  return false;
}
//****************************************************************************80

bool ppmb_write_header ( ofstream &output, int xsize, int ysize, int maxrgb )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_WRITE_HEADER writes the header of a binary portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &OUTPUT, a pointer to the file to contain the binary
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXRGB, the maximum RGB value.
//
//    Output, bool PPMB_WRITE_HEADER, is true if an error occurred.
//
{
  output << "P6"   << " " 
           << xsize  << " " 
           << ysize  << " " 
           << maxrgb << "\n";

  return false;
}
//****************************************************************************80

bool ppmb_write_test ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_WRITE_TEST tests the binary portable pixel map write routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to contain the binary
//    portable pixel map data.
//
//    Output, bool PPMB_WRITE_TEST is true if an error occurred.
//
{
  unsigned char *b;
  bool error;
  unsigned char *g;
  int maxrgb;
  unsigned char *r;
  int xsize;
  int ysize;

  xsize = 200;
  ysize = 200;
//
//  Allocate memory.
// 
  r = new unsigned char[ xsize * ysize ];
  g = new unsigned char[ xsize * ysize ];
  b = new unsigned char[ xsize * ysize ];
//
//  Set the data.
//
  error = ppmb_example ( xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_WRITE_TEST: Fatal error!\n";
    cout << "  PPM_EXAMPLE failed.\n";
    return true;
  }
//
//  Write the data to the file.
//
  error = ppmb_write ( file_name, xsize, ysize, r, g, b );

  delete [] r;
  delete [] g;
  delete [] b;

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_WRITE_TEST: Fatal error!\n";
    cout << "  PPMB_WRITE failed.\n";
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
