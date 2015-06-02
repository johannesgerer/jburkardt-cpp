# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "bmp_io.hpp"

int main ( long argc, char *argv[] );
void test01 ( char *file_name );
void test02 ( char *file_name );
void test03 ( char *file_name );
void test04 ( char *file_name );
void test05 ( char *file_name );
void timestamp ( );

//****************************************************************************80

int main ( long argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BMP_IO_PRB.
//
//  Discussion:
//
//    BMP_IO_PRB tests the BMP_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool bmp_byte_swap;
  char *file_name1 = "bmp_24.bmp";
  char *file_name2 = "blackbuck.bmp";
  char *file_name3 = "bmp_08.bmp";

  timestamp ( );
  cout << "\n";
  cout << "BMP_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BMP_IO library.\n";
//
//  Set the byte swapping option.
//
  bmp_byte_swap = true;
  bmp_byte_swap_set ( bmp_byte_swap );

  cout << "\n";
  cout << "  The BMP_BYTE_SWAP option being used is " << bmp_byte_swap << "\n";
//
//  Call BMP_WRITE_TEST to create a BMP file.
//
  test01 ( file_name1 );
  test02 ( file_name1 );
  test03 ( file_name1 );

  test04 ( file_name2 );

  test05 ( file_name3 );
//
//  Terminate.
//
  cout << "\n";
  cout << "BMP_IO_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 calls BMP_24_WRITE_TEST to create a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Call BMP_24_WRITE_TEST to store some graphics\n";
  cout << "  information in the file \"" << file_name << "\".\n";

  error = bmp_24_write_test ( file_name );

  if ( error )
  {
    cout << "\n";
    cout << "  BMP_24_WRITE_TEST failed.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  BMP_24_WRITE_TEST passed.\n";

  return;
}
//****************************************************************************80

void test02 ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 calls BMP_READ_TEST to read a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Call BMP_READ_TEST to read the\n";
  cout << "  information in the file \"" << file_name << "\".\n";

  error = bmp_read_test ( file_name );

  if ( error )
  {
    cout << "\n";
    cout << "  BMP_READ_TEST failed.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  BMP_READ_TEST passed.\n";

  return;
}
//****************************************************************************80

void test03 ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 calls BMP_PRINT_TEST to print a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Call BMP_PRINT_TEST to read and print the\n";
  cout << "  information in the file \"" << file_name << "\".\n";

  error = bmp_print_test ( file_name );

  if ( error )
  {
    cout << "\n";
    cout << "  BMP_PRINT_TEST failed.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  BMP_PRINT_TEST passed.\n";

  return;
}
//****************************************************************************80

void test04 ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 reads just the header from a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  unsigned char *barray;
  bool error;
  unsigned char *garray;
  long int height;
  unsigned char *rarray;
  unsigned long int width;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  BMP_READ can extract the RGB information from a BMP file.\n";
  cout << "\n";
  cout << "  We will try to read the file \"" << file_name << "\".\n";

  rarray = NULL;
  garray = NULL;
  barray = NULL;
//
//  Read the data from file.
//
  error = bmp_read ( file_name, &width, &height, &rarray, &garray,
    &barray );
//
//  Free the memory.
//
  delete [] rarray;
  delete [] garray;
  delete [] barray;

  if ( error )
  {
    cout << "\n";
    cout << "TEST04 - Fatal error!\n";
    cout << "  The test failed.\n";
  }
  else
  {
    cout << "\n";
    cout << "TEST04:\n";
    cout << "  WIDTH =  " << width  << ".\n";
    cout << "  HEIGHT = " << height << ".\n";
  }

  return;
}
//****************************************************************************80

void test05 ( char *file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 calls BMP_08_WRITE_TEST to create a BMP file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
{
  bool error;

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  Call BMP_08_WRITE_TEST to store some graphics\n";
  cout << "  information in the file \"" << file_name << "\".\n";

  error = bmp_08_write_test ( file_name );

  if ( error )
  {
    cout << "\n";
    cout << "  BMP_08_WRITE_TEST failed.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  BMP_08_WRITE_TEST passed.\n";

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
