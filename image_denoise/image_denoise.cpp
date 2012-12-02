# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "image_denoise.hpp"

//****************************************************************************80

int *gray_median_news ( int m, int n, int gray[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_MEDIAN_NEWS uses a median NEWS filter on a gray scale image to remove noise.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of pixels.
//
//    Input, int GRAY[M*N], the noisy grayscale data.
//
//    Output, int GRAY_MEDIAN_NEWS[M*N], the grayscale data for the filtered image.
//
{
  int *gray2;
  int i;
  int j;
  int p[5];

  gray2 = new int[m*n];
//
//  Process the main part of the image:
//
  for ( i = 1; i < m - 1; i++ )
  {
    for ( j = 1; j < n - 1; j++ )
    {
      p[0] = gray[i-1+ j   *m];
      p[1] = gray[i+1+ j   *m];
      p[2] = gray[i  +(j+1)*m];
      p[3] = gray[i  +(j-1)*m];
      p[4] = gray[i  + j   *m];

      gray2[i+j*m] = i4vec_median ( 5, p );
    }
  }
//
//  Process the four borders.
//  Get an odd number of data points, 
//
  for ( i = 1; i < m - 1; i++ )
  {
    j = 0;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i  + j   *m];
    p[3] = gray[i  +(j+1)*m];
    p[4] = gray[i  +(j+2)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );

    j = n - 1;
    p[0] = gray[i-1+ j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i  +(j-2)*m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  + j   *m];
    gray2[i+j*m] = i4vec_median ( 5, p );
  }

  for ( j = 1; j < n - 1; j++ )
  {
    i = 0;
    p[0] = gray[i  + j   *m];
    p[1] = gray[i+1+ j   *m];
    p[2] = gray[i+2+ j   *m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );

    i = m - 1;
    p[0] = gray[i-2+ j   *m];
    p[1] = gray[i-1+ j   *m];
    p[2] = gray[i  + j   *m];
    p[3] = gray[i  +(j-1)*m];
    p[4] = gray[i  +(j+1)*m];
    gray2[i+j*m] = i4vec_median ( 5, p );
  }
//
//  Process the four corners.
//
  i = 0;
  j = 0;
  p[0] = gray[i+1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j+1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = 0;
  j = n - 1;
  p[0] = gray[i+1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j-1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = m - 1;
  j = 0;
  p[0] = gray[i-1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j+1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  i = m - 1;
  j = n - 1;
  p[0] = gray[i-1+ j   *m];
  p[1] = gray[i  + j   *m];
  p[2] = gray[i  +(j-1)*m];
  gray2[i+j*m] = i4vec_median ( 3, p );

  return gray2;
}
//****************************************************************************80

int i4vec_frac ( int n, int a[], int k )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_FRAC searches for the K-th smallest entry in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    Hoare's algorithm is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, int A[N].
//    On input, A is the array to search.
//    On output, the elements of A have been somewhat rearranged.
//
//    Input, int K, the fractile to be sought.  If K = 1, the minimum
//    entry is sought.  If K = N, the maximum is sought.  Other values
//    of K search for the entry which is K-th in size.  K must be at
//    least 1, and no greater than N.
//
//    Output, double I4VEC_FRAC, the value of the K-th fractile of A.
//
{
  int frac;
  int i;
  int iryt;
  int j;
  int left;
  int temp;
  int x;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal nonpositive value of K = " << k << "\n";
    exit ( 1 );
  }

  if ( n < k )
  {
    cerr << "\n";
    cerr << "I4VEC_FRAC - Fatal error!\n";
    cerr << "  Illegal N < K, K = " << k << "\n";
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      frac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
//
//  Find I so that X <= A(I).
//
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
//
//  Find J so that A(J) <= X.
//
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp   = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return frac;
}
//****************************************************************************80

int i4vec_median ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MEDIAN returns the median of an unsorted I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    Hoare's algorithm is used.  The values of the vector are
//    rearranged by this routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input/output, int A[N], the array to search.  On output,
//    the order of the elements of A has been somewhat changed.
//
//    Output, int I4VEC_MEDIAN, the value of the median of A.
//
{
  int k;
  int median;

  k = ( n + 1 ) / 2;

  median = i4vec_frac ( n, a, k );

  return median;
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

void pgma_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxg )

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
//    Input, ifstream &FILE_IN, a pointer to the file.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXG, the maximum gray value.
//
{
  int count;
  char line[255];
  char *next;
  int step;
  int width;
  char word[255];

  step = 0;

  while ( 1 )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof() )
    {
      cerr << "\n";
      cerr << "PGMA_READ_HEADER - Fatal error!\n";
      cerr << "  End of file.\n";
      exit ( 1 );
    }

    next = line;

    if ( line[0] == '#' )
    {
      continue;
    }

    if ( step == 0 )
    {
      count = sscanf ( next, "%s%n", word, &width );
      if ( count == EOF )
      {
        continue;
      }
      next = next + width;
      if ( strcmp ( word, "P2" ) != 0 && strcmp ( word, "p2" ) != 0 )
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

      count = sscanf ( next, "%d%n", xsize, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      step = 2;
    }

    if ( step == 2 )
    {
      count = sscanf ( next, "%d%n", ysize, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      step = 3;
    }

    if ( step == 3 )
    {
      count = sscanf ( next, "%d%n", maxg, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      break;
    }

  }

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
