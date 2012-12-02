# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "bmp_io.hpp"
# include "ppma_io.hpp"

int main ( int argc, char *argv[] );
bool bmp_to_ppma ( char *filein_name, char *fileout_name );
void timestamp ( );
void ucmat_to_i4mat ( int xsize, int ysize, unsigned char *ucmat, int *i4mat );
void ucmat_vert_flip ( int xsize, int ysize, unsigned char *ucmat );

//****************************************************************************8080

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BMP_TO_PPMA.
//
//  Discussion:
//
//    BMP_TO_PPMA converts a BMP graphics file to an ASCII PPM file.
//
//    The program requires the BMP_IO and PPMA_IO libraries.
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
//    14 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    bmp_to_ppma file.bmp file.ppm
//
//  Parameters:
//
//    FILE.BMP is the name of the input BMP file to be read.
//
//    FILE.PPM is the name of the output ASCII PPM file to be created.
//
{
  bool error;
  char file_in_name[80];
  char file_out_name[80];

  timestamp ( );
  cout << "\n";
  cout << "BMP_TO_PPMA:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Convert a Microsoft BMP file to NetPBM ASCII PPM format.\n";
//
//  Get the specification for the input file.
//
  if ( argc < 2 )
  {
    cout << "\n";
    cout << "BMP_TO_PPMA:\n";
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
    cout << "BMP_TO_PPMA:\n";
    cout << "  Please enter the output PPMA file name:\n";

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
  error = bmp_to_ppma ( file_in_name, file_out_name );

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMA - Fatal error!\n";
    cout << "  The conversion was not successful.\n";
    exit ( 1 );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "BMP_TO_PPMA:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool bmp_to_ppma ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    BMP_TO_PPMA reads a BMP file and writes an ASCII PPM file.
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
//    Input, char *FILE_OUT_NAME, the name of the output ASCII PPM file.
//
//    Output, bool BMP_TO_PPMA, is true if an error occurred.
//
{
  unsigned char *b_ucmat;
  int *b_i4mat;
  bool error;
  unsigned char *g_ucmat;
  int *g_i4mat;
  long int height;
  unsigned char *r_ucmat;
  int *r_i4mat;
  int xsize;
  unsigned long int width;
  int ysize;

  r_ucmat = NULL;
  g_ucmat = NULL;
  b_ucmat = NULL;
//
//  Read the data from the BMP file.
//
  error = bmp_read ( file_in_name, &width, &height, &r_ucmat, &g_ucmat,
    &b_ucmat );

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMA - Fatal error!\n";
    cout << "  BMP_READ failed.\n";
    return error;
  }

  xsize = ( int ) width;
  ysize = ( int ) height;

  cout << "\n";
  cout << "BMP_TO_PPMA:\n";
  cout << "  XSIZE = " << xsize << ".\n";
  cout << "  YSIZE = " << ysize << ".\n";
//
//  The BMP up-down orientation is the opposite of what it is in PPM.
//
  cout << "\n";
  cout << "BMP_TO_PPMA:\n";
  cout << "  Flipping data orientation.\n";

  ucmat_vert_flip ( xsize, ysize, r_ucmat );
  ucmat_vert_flip ( xsize, ysize, g_ucmat );
  ucmat_vert_flip ( xsize, ysize, b_ucmat );
//
//  PPMA expects integers, not unsigned characters.
//
  r_i4mat = new int [xsize*ysize];
  g_i4mat = new int [xsize*ysize];
  b_i4mat = new int [xsize*ysize];

  ucmat_to_i4mat ( xsize, ysize, r_ucmat, r_i4mat );
  ucmat_to_i4mat ( xsize, ysize, g_ucmat, g_i4mat );
  ucmat_to_i4mat ( xsize, ysize, b_ucmat, b_i4mat );

  delete [] r_ucmat;
  delete [] g_ucmat;
  delete [] b_ucmat;
//
//  Write the PPMA file.
//
  error = ppma_write ( file_out_name, xsize, ysize,
    r_i4mat, g_i4mat, b_i4mat );

  if ( error )
  {
    cout << "\n";
    cout << "BMP_TO_PPMA - Fatal error!\n";
    cout << "  PPMA_WRITE failed.\n";
    return error;
  }

  delete [] r_i4mat;
  delete [] g_i4mat;
  delete [] b_i4mat;

  cout << "\n";
  cout << "BMP_TO_PPMA:\n";
  cout << "  Normal end of translation.\n";

  return false;
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

void ucmat_to_i4mat ( int xsize, int ysize, unsigned char *ucmat, int *i4mat )

//****************************************************************************80
//
//  Purpose:
//
//    UCMAT_TO_I4MAT converts a UCMAT to an I4MAT.
//
//  Discussion:
//
//    A UCMAT is a 2D array of unsigned characters.
//
//    An I4MAT is a 2D array of ints.
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
//    Input, unsigned char *UCMAT, the array of unsigned chars.
//
//    Output, int *I4MAT, the array of ints.
//
{
  int *i4mat2;
  unsigned char *ucmat2;
  int i;
  int j;
  unsigned char temp;

  ucmat2 = ucmat;
  i4mat2 = i4mat;

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      *i4mat2 = ( int ) *ucmat2;
      i4mat2 = i4mat2 + 1;
      ucmat2 = ucmat2 + 1;
    }
  }

  return;
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
