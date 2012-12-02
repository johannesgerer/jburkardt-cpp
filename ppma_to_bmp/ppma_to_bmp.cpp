# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "bmp_io.hpp"
# include "ppma_io.hpp"

int main ( int argc, char *argv[] );
void i4mat_vert_flip ( int xsize, int ysize, int *i4mat );
void i4mat_to_ucmat ( int xsize, int ysize, int *i4mat, unsigned char *ucmat );
bool ppma_to_bmp ( char *file_in_name, char *file_out_name );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PPMA_TO_BMP.
//
//  Discussion:
//
//    PPMA_TO_BMP converts an ASCII PPM file to a BMP file.
//
//    The program requires the BMP_IO and PPMA_IO libraries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    ppma_to_bmp file.ppma file.bmp
//
//  Parameters:
//
//    FILE.PPMA is the name of the input ASCII PPM file to be read.
//
//    FILE.BMP is the name of the output BMP file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );

  cout << "\n";
  cout << "PPMA_TO_BMP:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a NetPBM ASCII PPM file to Microsoft BMP format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PPMA_TO_BMP:\n";
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
    cout << "PPMA_TO_BMP:\n";
    cout << "  Please enter the output BMP file name:\n";

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

  error = ppma_to_bmp ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_BMP - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "PPMA_TO_BMP:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool ppma_to_bmp ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_TO_BMP reads a PPMA (ASCII PPM) file and writes a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input PPMA file.
//
//    Input, char *FILE_OUT_NAME, the name of the output BMP file.
//
//    Output, bool PPMA_TO_BMP, is true if an error occurred.
//
{
  int *b_i4mat;
  unsigned char *b_ucmat;
  bool error;
  int *g_i4mat;
  unsigned char *g_ucmat;
  long int height;
  int maxrgb;
  int *r_i4mat;
  unsigned char *r_ucmat;
  unsigned long int width;
  int xsize;
  int ysize;
//
//  Read the integer data from the PPMA file.
//
  error = ppma_read ( file_in_name, xsize, ysize, maxrgb,
    &r_i4mat, &g_i4mat, &b_i4mat );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_BMP - Fatal error!\n";
    cout << "  PPMA_READ failed.\n";
    return error;
  }

  cout << "\n";
  cout << "PPMA_TO_BMP:\n";
  cout << "  XSIZE = " << xsize << ".\n";
  cout << "  YSIZE = " << ysize << ".\n";
//
//  The BMP up-down orientation is the opposite of what it is in PPM.
//
  cout << "\n";
  cout << "PPMA_TO_BMP:\n";
  cout << "  Flipping data orientation.\n";

  i4mat_vert_flip ( xsize, ysize, r_i4mat );
  i4mat_vert_flip ( xsize, ysize, g_i4mat );
  i4mat_vert_flip ( xsize, ysize, b_i4mat );
//
//  BMP expects unsigned characters, not integers.
//
  r_ucmat = new unsigned char [xsize*ysize];
  g_ucmat = new unsigned char [xsize*ysize];
  b_ucmat = new unsigned char [xsize*ysize];

  i4mat_to_ucmat ( xsize, ysize, r_i4mat, r_ucmat );
  i4mat_to_ucmat ( xsize, ysize, g_i4mat, g_ucmat );
  i4mat_to_ucmat ( xsize, ysize, b_i4mat, b_ucmat );

  delete [] r_i4mat;
  delete [] g_i4mat;
  delete [] b_i4mat;

  height = ( long int ) ysize;
  width = ( unsigned long int ) xsize;

  error = bmp_24_write ( file_out_name, width, height,
    r_ucmat, g_ucmat, b_ucmat );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_TO_BMP - Fatal error!\n";
    cout << "  BMP_08_WRITE failed.\n";
    return error;
  }

  delete [] r_ucmat;
  delete [] g_ucmat;
  delete [] b_ucmat;

  cout << "\n";
  cout << "PPMA_TO_BMP:\n";
  cout << "  Normal end of translation.\n";

  return false;
}
//****************************************************************************80

void i4mat_vert_flip ( int xsize, int ysize, int *i4mat )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_VERT_FLIP swaps rows of an I4MAT, to flip it vertically.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, int *I4MAT, the address of the first element of the array.
//
{
  int i;
  int j;
  int k1;
  int k2;
  int temp;

  for ( j = 0; j <= ( ysize / 2 ); j++ )
  {
    k1 = xsize * j;
    k2 = xsize * ( ysize - 1 - j );

    for ( i = 0; i < xsize; i++ )
    {
      temp        = i4mat[k1+i];
      i4mat[k1+i] = i4mat[k2+i];
      i4mat[k2+i] = temp;
    }
  }

  return;
}
//****************************************************************************80

void i4mat_to_ucmat ( int xsize, int ysize, int *i4mat, unsigned char *ucmat )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TO_UCMAT converts an I4MAT to a UCMAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, int char *I4MAT, the array of ints.
//
//    Output, unsigned char *UCMAT, the array of unsigned chars.
//
{
  int *i4mat2;
  unsigned char *ucmat2;
  int i;
  int j;

  ucmat2 = ucmat;
  i4mat2 = i4mat;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      *ucmat2 = ( unsigned char ) *i4mat;
      i4mat2 = i4mat2 + 1;
      ucmat2 = ucmat2 + 1;
    }
  }

  return;
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
