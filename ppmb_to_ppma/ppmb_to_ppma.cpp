# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
int i4_max ( int i1, int i2 );
bool ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
  int *g, int *b );
bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *r,
  int *g, int *b );
bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize,
  int ysize, int rgb_max );
bool ppmb_check_data ( int xsize, int ysize, int maxrgb, unsigned char *r,
  unsigned char *g, unsigned char *b );
bool ppmb_read ( char *file_in_name, int *xsize, int *ysize, int *maxrgb,
  unsigned char **r, unsigned char **g, unsigned char **b );
bool ppmb_read_data ( ifstream &file_in, int xsize, int ysize,
  unsigned char *r, unsigned char *g, unsigned char *b );
bool ppmb_read_header ( ifstream &file_in, int *xsize, int *ysize,
  int *maxrgb );
bool ppmb_to_ppma ( char *file_in_name, char *file_out_name );
void timestamp ( );
void ucmat_to_i4mat ( int n, unsigned char *a, int *b );

//****************************************************************************80*

int main ( int argc, char *argv[] )

//****************************************************************************80*
//
//  Purpose:
//
//    MAIN is the main program for PPMB_TO_PPMA.
//
//  Discussion:
//
//    PPMB_TO_PPMA converts a binary PPM file to an ASCII PPM file.
//
//    The application requires the PPMA_IO and PPMB_IO libraries.
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
//  Usage:
//
//    ppmb_to_ppma file.ppmb file.ppma
//
//  Parameters:
//
//    FILE.PPMB is the name of the input binary PPM file to be read.
//
//    FILE.PPMA is the name of the output ASCII PPM file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );
  cout << "\n";
  cout << "PPMB_TO_PPMA:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a NetPBM binary PPM file to ASCII PPM format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PPMB_TO_PPMA:\n";
    cout << "  Please enter the input PPMB file name:\n";

    cin.getline ( file_in_name, sizeof ( file_in_name ) );

    error = cin.eof ( );

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
    cout << "PPMB_TO_PPMA:\n";
    cout << "  Please enter the output PPMA file name:\n";

    cin.getline ( file_out_name, sizeof ( file_out_name ) );

    error = cin.eof ( );

    if ( error )
    {
      exit ( 1 );
    }
  }
  else
  {
    strcpy ( file_out_name, argv[2] );
  }

  error = ppmb_to_ppma ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_PPMA - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "PPMB_TO_PPMA:\n";
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

bool ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
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
//    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
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
  ofstream file_out;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
//
//  Open the output file.
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\".\n";
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
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, rgb_max );

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
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

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
  file_out.close ( );

  return false;
}
//****************************************************************************80

bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *r,
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
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
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
      file_out << *r_index << " " << *g_index << " " << *b_index;
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        file_out << "\n";
      }
      else
      {
        file_out << " ";
      }
    }
  }
  return false;
}
//****************************************************************************80

bool ppma_write_header ( ofstream &file_out, char *file_out_name, int xsize,
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
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, char *FILE_OUT_NAME, the name of the file.
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
  file_out << "P3\n";
  file_out << "# " << file_out_name << " created by PPMA_IO::PPMA_WRITE.C.\n";
  file_out << xsize << "  " << ysize << "\n";
  file_out << rgb_max << "\n";

  return false;
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

bool ppmb_read ( char *file_in_name, int *xsize, int *ysize, int *maxrgb,
  unsigned char **r, unsigned char **g, unsigned char **b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_READ reads the header and data from a binary portable pixel map file.
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
//    Input, char *FILE_IN_NAME, the name of the file containing the binary
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXRGB, the maximum RGB value.
//
//    Output, unsigned char **R, **G, **B, the arrays of XSIZE by YSIZE
//    data values.
//
//    Output, bool PPMB_READ, is true if there was an error.
//
{
  bool error;
  ifstream file_in;

  file_in.open ( file_in_name, ios::binary );

  if ( !file_in )
  {
    cout << "\n";
    cout << "PPMB_READ - Fatal error!\n";
    cout << "  Cannot open the input file \"" << file_in_name << "\".\n";
    return true;
  }
//
//  Read the header.
//
  error = ppmb_read_header ( file_in, xsize, ysize, maxrgb );

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
  *r = new unsigned char[(*xsize) * (*ysize)];
  *g = new unsigned char[(*xsize) * (*ysize)];
  *b = new unsigned char[(*xsize) * (*ysize)];
//
//  Read the data.
//
  error = ppmb_read_data ( file_in, *xsize, *ysize, *r, *g, *b );

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
  file_in.close ( );

  return false;
}
//****************************************************************************80

bool ppmb_read_data ( ifstream &file_in, int xsize, int ysize,
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
//    Input, ifstream &file_in, a pointer to the file containing the binary
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
      file_in.read ( &c, 1 );
      *indexr = ( unsigned char ) c;
      indexr = indexr + 1;
      error = file_in.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PPMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading R byte of pixel ("
          << i << ", " << j <<") \n";
        return true;
      }

      file_in.read ( &c, 1 );
      *indexg = ( unsigned char ) c;
      indexg = indexg + 1;
      error = file_in.eof();
      if ( error )
      {
        cout << "\n";
        cout << "PPMB_READ_DATA - Fatal error!\n";
        cout << "  End of file reading G byte of pixel ("
          << i << ", " << j <<") \n";
        return true;
      }

      file_in.read ( &c, 1 );
      *indexb = ( unsigned char ) c;
      indexb = indexb + 1;
      error = file_in.eof();
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

bool ppmb_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxrgb )

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
//    12 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &file_in, a pointer to the file containing the binary
//    portable pixel map data.
//
//    Output, int *XSIZE, *YSIZE, the number of rows and columns of data.
//
//    Output, int *MAXRGB, the maximum RGB value.
//
//    Output, bool PPMB_READ_HEADER, is true if an error occurred.
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
      cout << "PPMB_READ_HEADER - Fatal error!\n";
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
      if ( strcmp ( word, "P6" ) != 0 && strcmp ( word, "p6" ) != 0 )
      {
        cout << "\n";
        cout << "PPMB_READ_HEADER - Fatal error.\n";
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
      count = sscanf ( next, "%d%n", maxrgb, &width );
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
//****************************************************************************80*/

bool ppmb_to_ppma ( char *file_in_name, char *file_out_name )

//****************************************************************************80*/
//
//  Purpose:
//
//    PPMB_TO_PPMA converts one PPMB file to PPMA format.
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
//    Input, char *FILE_IN_NAME, the name of the input PPMB file.
//
//    Input, char *FILE_OUT_NAME, the name of the output PPMA file.
//
//    Output, bool PPMB_TO_PPMA, is true if an error occurred.
//
{
  unsigned char* b;
  int *b2;
  bool error;
  unsigned char* g;
  int *g2;
  int maxrgb;
  unsigned char* r;
  int *r2;
  int xsize;
  int ysize;
//
//  Read the input file.
//
  error = ppmb_read ( file_in_name, &xsize, &ysize, &maxrgb, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_PPMA: Fatal error!\n";
    cout << "  PPMB_READ failed.\n";
    return true;
  }
//
//  Check the data.
//
  error = ppmb_check_data ( xsize, ysize, maxrgb, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_PPMA: Fatal error!\n";
    cout << "  PPMB_CHECK_DATA reports bad data from the file.\n";

    delete [] r;
    delete [] g;
    delete [] b;
    return true;
  }
//
//  Convert the data.
//
  r2 = new int [ xsize *  ysize ];
  ucmat_to_i4mat ( xsize * ysize, r, r2 );
  delete [] r;

  g2 = new int [ xsize *  ysize ];
  ucmat_to_i4mat ( xsize * ysize, g, g2 );
  delete [] g;

  b2 = new int [ xsize *  ysize ];
  ucmat_to_i4mat ( xsize * ysize, b, b2 );
  delete [] b;
//
//  Write the output file.
//
  error = ppma_write ( file_out_name, xsize, ysize, r2, g2, b2 );

  delete [] r2;
  delete [] g2;
  delete [] b2;

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_PPMA: Fatal error!\n";
    cout << "  PPMA_WRITE failed.\n";
    return true;
  }

  cout << "\n";
  cout << "PPMB_TO_PPMA:\n";
  cout << "  Normal end of execution.\n";

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
//****************************************************************************80

void ucmat_to_i4mat ( int n, unsigned char *a, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    UCMAT_TO_I4MAT converts a UCMAT to an I4MAT.
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
//    Input, unsigned int *B, a pointer to a vector of ints.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    *b = ( int ) *a;
    a++;
    b++;
  }

}
