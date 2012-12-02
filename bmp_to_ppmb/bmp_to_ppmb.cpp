# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "bmp_io.hpp"
# include "ppmb_io.hpp"

int main ( int argc, char *argv[] );
int bmp_to_ppmb ( char *filein_name, char *fileout_name );
void timestamp ( );
void ucmat_vert_flip ( int xsize, int ysize, unsigned char *ucmat );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BMP_TO_PPMB.
//
//  Discussion:
//
//    BMP_TO_PPMB converts a BMP graphics file to a binary PPM file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Modified:
//
//    14 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    bmp_to_ppmb file.bmp file.ppm
//
//  Parameters:
//
//    FILE.BMP is the name of the input BMP file to be read.
//
//    FILE.PPM is the name of the output binary PPM file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );
  cout << "\n";
  cout << "BMP_TO_PPMB:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a Microsoft BMP file to NetPBM binary PPM format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "BMP_TO_PPMB:\n";
    cout << "  Please enter the input BMP file name:\n";

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
    cout << "BMP_TO_PPMB:\n";
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

  error = bmp_to_ppmb ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMB - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "BMP_TO_PPMB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int bmp_to_ppmb ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    BMP_TO_PPMB reads a BMP file and writes a PPMB file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Modified:
//
//    13 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input BMP file.
//
//    Input, char *FILE_OUT_NAME, the name of the output PPMA file.
//
{
  unsigned char *b;
  bool error;
  unsigned char *g;
  long int height;
  unsigned char *r;
  unsigned long int width;
  int xsize;
  int ysize;
//
//  Read the data from the BMP file.
//
  error = bmp_read ( file_in_name, &width, &height, &r, &g, &b );

  xsize = ( int ) width;
  ysize = ( int ) height;

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMB - Fatal error!\n";
    cout << "  BMP_READ failed.\n";
    return error;
  }
//
//  The BMP up-down orientation is the opposite of what it is in PPM.
//
  ucmat_vert_flip ( xsize, ysize, r );
  ucmat_vert_flip ( xsize, ysize, g );
  ucmat_vert_flip ( xsize, ysize, b );
//
//  Write the PPMB file.
//
  error = ppmb_write ( file_out_name, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMB - Fatal error!\n";
    cout << "  PPMB_WRITE failed.\n";
    return error;
  }

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
//    29 April 2003
//
//  Author:
//
//    John Burkardt
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
//  Discussion:
//
//    A UCMAT is a 2D array of unsigned characters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int XSIZE, int YSIZE, the number of columns and rows.
//
//    Input, unsigned char *UCMAT, the address of the first element
//    of the array.
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
