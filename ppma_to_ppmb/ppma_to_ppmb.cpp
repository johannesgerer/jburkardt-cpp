# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int i4_max ( int i1, int i2 );
void i4vec_to_ucvec ( int n, int *i4vec, unsigned char *ucvec );
bool ppma_check_data ( int xsize, int ysize, int rgb_max, int *r,
  int *g, int *b );
bool ppma_read ( char *file_in_name, int *xsize, int *ysize, int *rgb_max,
  int **r, int **g, int **b );
bool ppma_read_data ( ifstream &file_in, int xsize, int ysize, int *r,
  int *g, int *b );
bool ppma_read_header ( ifstream &file_in, int *xsize, int *ysize,
  int *rgb_max );
bool ppma_to_ppmb ( char *file_in_name, char *file_out_name );
bool ppmb_write ( char *file_out_name, int xsize, int ysize,
  unsigned char *r, unsigned char *g, unsigned char *b );
bool ppmb_write_data ( ofstream &file_out, int xsize, int ysize,
  unsigned char *r, unsigned char *g, unsigned char *b );
bool ppmb_write_header ( ofstream &file_out, int xsize, int ysize, int maxrgb );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PPMA_TO_PPMB.
//
//  Discussion:
//
//    PPMA_TO_PPMB converts an ASCII PPM file to a binary PPM file.
//
//    This program requires the PPMA_IO and PPMB_IO libraries.
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
//  Usage:
//
//    ppma_to_ppmb file.ppma file.ppmb
//
//  Parameters:
//
//    FILE.PPMA is the name of the input ASCII PPM file to be read.
//
//    FILE.PPMB is the name of the output binary PPM file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );
  cout << "\n";
  cout << "PPMA_TO_PPMB:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a NetPBM ASCII PPM file to binary PPM format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PPMA_TO_PPMB:\n";
    cout << "  Please enter the input PPMA file name:\n";

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
    cout << "PPMA_TO_PPMB:\n";
    cout << "  Please enter the output PPMB file name:\n";

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

  error = ppma_to_ppmb ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_PPMB - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "PPMA_TO_PPMB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
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

void i4vec_to_ucvec ( int n, int *i4vec, unsigned char *ucvec )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_TO_UCVEC converts an I4VEC to a UCVEC.
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
//    Input, int *I4VEC, a pointer to a vector of ints.
//
//    Input, unsigned char *UCVEC, a pointer to a vector of unsigned chars.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    *ucvec = ( unsigned char ) *i4vec;
    i4vec++;
    ucvec++;
  }

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
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  xsize <= 0.\n";
    cout << "  xsize = " << xsize << "\n";
    return true;
  }

  if ( ysize <= 0 )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  ysize <= 0.\n";
    cout << "  ysize = " << ysize << "\n";
    return true;
  }

  if ( r == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  Null pointer to R.\n";
    return true;
  }

  if ( g == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
    cout << "  Null pointer to G.\n";
    return true;
  }

  if ( b == NULL )
  {
    cout << "\n";
    cout << "PPMA_CHECK_DATA: Error!\n";
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

bool ppma_read ( char *file_in_name, int *xsize, int *ysize, int *rgb_max,
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
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the file containing the ASCII
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *RGB_MAX, the maximum RGB value.
//
//    Output, int **R, **G, **B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_READ, is
//    true, if an error was detected, or
//    false, if the file was read.
//
{
  bool error;
  ifstream file_in;
  int numbytes;

  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "PPMA_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << file_in_name << "\".\n";
    return true;
  }
//
//  Read the header.
//
  error = ppma_read_header ( file_in, xsize, ysize, rgb_max );

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
  numbytes = (*xsize) * (*ysize) * sizeof ( int );

  *r = new int[numbytes];
  *g = new int[numbytes];
  *b = new int[numbytes];
//
//  Read the data.
//
  error = ppma_read_data ( file_in, *xsize, *ysize, *r, *g, *b );

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
  file_in.close ( );

  return false;
}
//****************************************************************************80

bool ppma_read_data ( ifstream &file_in, int xsize, int ysize, int *r,
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
//    Input, ifstream &FILE_IN, a pointer to the file containing the ASCII
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
      file_in >> *r;
      if ( file_in.eof() )
      {
        return true;
      }
      r = r + 1;

      file_in >> *g;
      if ( file_in.eof() )
      {
        return true;
      }
      g = g + 1;

      file_in >> *b;
      if ( file_in.eof() )
      {
        return true;
      }
      b = b + 1;
    }
  }

  return false;
}
//****************************************************************************80

bool ppma_read_header ( ifstream &file_in, int *xsize, int *ysize,
  int *rgb_max )

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
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_IN, a pointer to the file containing the ASCII
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *RGB_MAX, the maximum RGB value.
//
//    Output, bool PPMA_READ_HEADER, is
//    true, if an error was detected, or
//    false, if the header was read.
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
      cout << "\n";
      cout << "PPMA_READ_HEADER - Fatal error!\n";
      cout << "  End of file.\n";
      return true;
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
      if ( strcmp ( word, "P3" ) != 0 && strcmp ( word, "p3" ) != 0 )
      {
        cout << "\n";
        cout << "PPMA_READ_HEADER - Fatal error.\n";
        cout << "  Bad magic number = \"" << word << "\".\n";
        return true;
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
      count = sscanf ( next, "%d%n", rgb_max, &width );
      next = next + width;
      if ( count == EOF )
      {
        continue;
      }
      break;
    }

  }

  return false;
}
//****************************************************************************80

bool ppma_to_ppmb ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_TO_PPMB converts an ASCII PPM file to binary PPM format.
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
//    Input, char *FILE_IN_NAME, the name of the input ASCII PPM file to be read.
//
//    Input, char *FILE_OUT_NAME, the name of the output binary PPM file to be created.
//
//    Output, bool PPMA_TO_PPMB, is true if an error occurred.
//
{
  int *b;
  unsigned char *b2;
  int *g;
  unsigned char *g2;
  bool error;
  int  maxrgb;
  int *r;
  unsigned char *r2;
  int  xsize;
  int  ysize;
//
//  Read the input file.
//
  error = ppma_read ( file_in_name, &xsize, &ysize, &maxrgb, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_PPMB: Fatal error!\n";
    cout << "  PPMA_READ failed.\n";
    return true;
  }
//
//  Check the data.
//
  error = ppma_check_data ( xsize, ysize, maxrgb, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_PPMB: Fatal error!\n";
    cout << "  PPM_CHECK_DATA reports bad data from the file.\n";
    return true;
  }
//
//  Copy the data into unsigned char's.
//
  r2 = new unsigned char [ xsize * ysize ];
  i4vec_to_ucvec ( xsize * ysize, r, r2 );
  delete [] r;

  g2 = new unsigned char [ xsize * ysize ];
  i4vec_to_ucvec ( xsize * ysize, g, g2 );
  delete [] g;

  b2 = new unsigned char [ xsize * ysize ];
  i4vec_to_ucvec ( xsize * ysize, b, b2 );
  delete [] b;
//
//  Write the output file.
//
  error = ppmb_write ( file_out_name, xsize, ysize, r2, g2, b2 );

  delete [] r2;
  delete [] g2;
  delete [] b2;

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_PPMB: Fatal error!\n";
    cout << "  PPMB_WRITE failed.\n";
    return true;
  }

  return false;
}
//****************************************************************************80

bool ppmb_write ( char *file_out_name, int xsize, int ysize,
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
//    Input, char *FILE_OUT_NAME, the name of the file to contain the binary
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
  ofstream file_out;
  int i;
  unsigned char *indexb;
  unsigned char *indexg;
  unsigned char *indexr;
  int j;
  int maxrgb;
//
//  Open the output file.
//
  file_out.open ( file_out_name, ios::binary );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PPMB_WRITE: Fatal error!\n";
    cout << "  Cannot open the output file " << file_out_name << ".\n";
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
  error = ppmb_write_header ( file_out, xsize, ysize, maxrgb );

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
  error = ppmb_write_data ( file_out, xsize, ysize, r, g, b );

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
  file_out.close ( );

  return false;
}
//****************************************************************************80

bool ppmb_write_data ( ofstream &file_out, int xsize, int ysize,
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
//    Input, ofstream FILE_OUT, a pointer to the file to contain the binary
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
      file_out << *indexr;
      indexr = indexr + 1;
      file_out << *indexg;
      indexg = indexg + 1;
      file_out << *indexb;
      indexb = indexb + 1;
    }
  }
  return false;
}
//****************************************************************************80

bool ppmb_write_header ( ofstream &file_out, int xsize, int ysize, int maxrgb )

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
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the binary
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int MAXRGB, the maximum RGB value.
//
//    Output, bool PPMB_WRITE_HEADER, is true if an error occurred.
//
{
  file_out << "P6"   << " "
           << xsize  << " "
           << ysize  << " "
           << maxrgb << "\n";

  return false;
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
