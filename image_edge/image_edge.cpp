# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "image_edge.hpp"

//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
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
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
}
//****************************************************************************80

int *i4mat_histogram ( int m, int n, int a[], int histo_num )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//    It is assumed that the entries in the vector A are nonnegative.
//    Only values between 0 and HISTO_NUM will be histogrammed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, int A[M*N], the array to examine.
//
//    Input, int HISTO_NUM, the maximum value for which a
//    histogram entry will be computed.
//
//    Output, int I4MAT_HISTOGRAM[HISTO_NUM+1], contains the number of
//    entries of A with the values of 0 through HISTO_NUM.
//
{
  int *histo_gram;
  int i;
  int j;

  histo_gram = new int[histo_num+1];

  for ( i = 0; i <= histo_num; i++ )
  {
    histo_gram[i] = 0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( 0 <= a[i+j*m] && a[i+j*m] <= histo_num )
      {
        histo_gram[a[i+j*m]] = histo_gram[a[i+j*m]] + 1;
      }
    }
  }

  return histo_gram;
}
//****************************************************************************80

int i4mat_max ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MAX returns the maximum of an I4MAT.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
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
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Output, int I4MAT_MAX, the maximum entry of A.
//
{
  int i;
  int j;
  int value;

  value = - i4_huge ( );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

int *news ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    NEWS demonstrates the NEWS stencil for edge detection.
//
//  Discussion:
//
//    Given a black and white image A, which we regard as an M by N array
//    of pixels, we want to produce an array E of the same shape, which
//    contains information describing the location of edges.
//
//    A simple algorithm for trying to detect edges in an array that
//    represents an image is the NEWS scheme.  For each pixel A(C),
//    we consider its North, East, West, and South pixel neighbors.  The
//    indexing of arrays and images do not correspond, so we will use
//    these directions instead:
//
//             A(N)
//              |
//              |
//      A(W)---A(C)---A(E)
//              |
//              |
//             A(S)
//
//    Entry E(C) of the edge array will be computed by
//
//      E(C) = abs ( A(N) - A(S) ) + abs ( A(E) - A(W) )
//
//    Pixels of A that represent edges will tend to have high values
//    of E, while pixels that are interior to a region of roughly the
//    same shade will tend to have low values.
//
//    Thus, an edge detection scheme would use the NEWS stencil to
//    compute the E array, determine E_MAX, the maximum entry in E,
//    choose some threshold value E_THRESH, and declare pixel A(I,J)
//    to be associated with an edge whenever E(I,J) is greater than E_THRESH.
//
//    In this program, we demonstrate the NEWS stencil using a PGM
//    grayscale image of coins.  At the end, we use the edge information
//    to produce a color image in which the edges of the coins have been
//    outlined in red.
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
//    Input, int M, N, the number of rows and columns in the image.
//
//    Input, int A[M*N], the gray scale image data, presumably
//    integers between 0 and 255.
//
//    Output, int NEWS[M*N], is 1 for each pixel that is part of an
//    edge, and 0 otherwise.
//
{
  int *b;
  int *e;
  int e_max;
  int i;
  int j;
  int thresh;
//
//  For neatness, we add a border of zeros to the image,
//  then fill in the border by copying the nearby original values.
//  This will be our M+2 by N+2 data array B.
//
  b = new int[(m+2)*(n+2)];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[i+1+(j+1)*(m+2)] = a[i+j*m];
    }
  }

  for ( j = 1; j < n + 1; j++ )
  {
    b[  0+j*(m+2)] = b[1+j*(m+2)];
    b[m+1+j*(m+2)] = b[m+j*(m+2)];
  }

  for ( i = 1; i < m + 1; i++ )
  {
    b[i+   0 *(m+2)] = b[i+1*(m+2)];
    b[i+(n+1)*(m+2)] = b[i+n*(m+2)];
  }

  b[  0+   0 *(m+2)] = ( b[  0+   1 *(m+2)] + b[  1+   0 *(m+2)] ) / 2;
  b[m+1+   0 *(m+2)] = ( b[m+1+   1 *(m+2)] + b[m+0+   0 *(m+2)] ) / 2;
  b[  0+(n+1)*(m+2)] = ( b[  0+(n+0)*(m+2)] + b[  1+(n+1)*(m+2)] ) / 2;
  b[m+1+(n+1)*(m+2)] = ( b[m+1+(n+0)*(m+2)] + b[m+0+(n+1)*(m+2)] ) / 2;
//
//  Apply the NEWS Operator.  We do not process the boundary pixels.
//
//  The picture is:
//
//   |  0 +1  0 |     |  0  0   0 |
//   |  0  0  0 |  +  | -1  0  +1 |
//   |  0 -1  0 |     |  0  0   0 |
//
  e = new int[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      e[i+j*m] = abs ( - b[i  +(j+1)*(m+2)] + b[i+2+(j+1)*(m+2)] )
               + abs ( - b[i+1+ j   *(m+2)] + b[i+1+(j+2)*(m+2)] );
    }
  }
  delete [] b;
//
//  Remap E so the largest value is 255.
//
  e_max = i4mat_max ( m, n, e );
//
//  Threshold the data.  Set the threshold to give enough detail
//  to guess the coin denominations.
//
  thresh = e_max / 5;

  cout << "\n";
  cout << "NEWS:\n";
  cout << "  E_MAX = " << e_max << "\n";
  cout << "  Using threshold value THRESH = " << thresh << "\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( e[i+j*m] < thresh )
      {
        e[i+j*m] = 0;
      }
      else
      {
        e[i+j*m] = 1;
      }
    }
  }
  return e;
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

void pgma_read_data ( ifstream &file_in, int xsize, int ysize, int *g )

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
//    Input, ifstream &FILE_IN, a pointer to the file containing the data.
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
      file_in >> *g;
      if ( file_in.eof ( ) )
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
