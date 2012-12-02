# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "bmp_io.hpp"
# include "ppmb_io.hpp"

int main ( int argc, char *argv[] );
bool ppmb_to_bmp ( char *file_in_name, char *file_out_name );
void timestamp ( );
void ucmat_vert_flip ( int xsize, int ysize, unsigned char *ucmat );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PPMB_TO_BMP.
//
//  Discussion:
//
//    PPMB_TO_BMP converts a binary PPM file to a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    ppmb_2_bmp file.ppmb file.bmp
//
//  Parameters:
//
//    FILE.PPMB is the name of the input binary PPM file to be created.
//
//    FILE.BMP is the name of the output BMP file to be read.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );

  cout << "\n";
  cout << "PPMB_TO_BMP:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a NetPBM binary PPM file to Microsoft BMP format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "PPMB_TO_BMP:\n";
    cout << "  Please enter the input PPMB file name:\n";

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
    cout << "PPMB_TO_BMP:\n";
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

  error = ppmb_to_bmp ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_BMP - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PPMB_TO_BMP:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool ppmb_to_bmp ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PPMB_TO_BMP reads a PPMB (binary PPM) file and writes a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 December 2011
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
{
  unsigned char *b;
  bool error;
  unsigned char *g;
  long int height;
  int maxrgb;
  unsigned char *r;
  unsigned long int width;
  int xsize;
  int ysize;
//
//  Read the data from the PPMB file.
//
  error = ppmb_read ( file_in_name, xsize, ysize, maxrgb, &r, &g, &b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_BMP - Fatal error!\n";
    cout << "  PPMB_READ failed.\n";
    return true;
  }

  cout << "\n";
  cout << "PPMB_TO_BMP:\n";
  cout << "  XSIZE = " << xsize << "\n";
  cout << "  YSIZE = " << ysize << "\n";
//
//  The BMP up-down orientation is the opposite of what it is in PPM.
//
  ucmat_vert_flip ( xsize, ysize, r );
  ucmat_vert_flip ( xsize, ysize, g );
  ucmat_vert_flip ( xsize, ysize, b );

  width = ( unsigned long int ) xsize;
  height = ( long int ) ysize;

  error = bmp_24_write ( file_out_name, width, height, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMB_TO_BMP - Fatal error!\n";
    cout << "  BMP_WRITE failed.\n";
    return true;
  }
//
//  Free the memory.
//
  delete [] r;
  delete [] g;
  delete [] b;

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

void ucmat_vert_flip ( int xsize, int ysize, unsigned char *ucmat )

//****************************************************************************80
//
//  Purpose:
//
//    UCMAT_VERT_FLIP swaps rows of a UCMAT, to flip it vertically.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, unsigned char *UCMAT, the address of the first element of the array.
//
{
  int i;
  int j;
  int k1;
  int k2;
  unsigned char temp;

  for ( j = 0; j <= ( ysize / 2 ); j++ )
  {
    k1 = xsize * j;
    k2 = xsize * ( ysize - 1 - j );

    for ( i = 0; i < xsize; i++ )
    {
      temp        = ucmat[k1+i];
      ucmat[k1+i] = ucmat[k2+i];
      ucmat[k2+i] = temp;
    }
  }

  return;
}
