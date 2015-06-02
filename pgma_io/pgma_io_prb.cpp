# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "pgma_io.hpp"

int main ( int argc, char *argv[] );
void test01 ( );
void test02 ( );
void test03 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PGMA_IO_PRB.
//
//  Discussion:
//
//    PGMA_IO_PRB tests the PGMA_IO library.
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
{
  bool error;

  timestamp ( );
  cout << "\n";
  cout << "PGMA_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PGMA_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PGMA_IO_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )  

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PGMA_EXAMPLE, PGMA_WRITE.
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
{
  string output_name = "pgma_io_prb_01.ascii.pgm";
  int *g;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  PGMA_EXAMPLE sets up ASCII PGM data.\n";
  cout << "  PGMA_WRITE writes an ASCII PGM file.\n";
  cout << "\n";
  cout << "  Writing the file \"" << output_name << "\".\n";

  g = new int[xsize*ysize];

  pgma_example ( xsize, ysize, g );

  cout << "\n";
  cout << "  PGMA_EXAMPLE has set up the data.\n";

  pgma_write ( output_name, xsize, ysize, g );

  cout << "\n";
  cout <<  "  PGMA_WRITE was successful.\n";

  delete [] g;
//
//  Now have PGMA_READ_TEST look at the file we think we created.
//
  pgma_read_test ( output_name );

  cout << "\n";
  cout <<  "  PGMA_READ_TEST was able to read our file.\n";

  return;
}
//****************************************************************************80

void test02 ( )  

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PGMA_READ_HEADER, PGMA_READ_DATA.
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
{
  string input_name = "pgma_io_prb_02.ascii.pgm";
  int *g;
  int i;
  int j;
  int k;
  int maxg;
  int xsize;
  int ysize;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PGMA_READ reads the header and data of an ASCII PGM file.\n";
  cout << "\n";
  cout << "  Reading the file \"" << input_name << "\".\n";
//
//  Create a data file to read.
//
  pgma_write_test ( input_name );

  cout << "\n";
  cout << "  PGMA_WRITE_TEST created some test data.\n";
//
//  Now have PGMA_READ try to read it.
//
  pgma_read ( input_name, xsize, ysize, maxg, &g );

  cout << "\n";
  cout << "  PGMA_READ read the test data successfully.\n";

  cout << "\n";
  cout << "  Sample data:\n";
  cout << "\n";
  for ( k = 0; k <= 9; k++ )
  {
    i = ( ( 9 - k ) * 0 + k * ( xsize - 1 ) ) / 9;
    j = ( ( 9 - k ) * 0 + k * ( ysize - 1 ) ) / 9;
    cout << setw(4) << i            << "  "
         << setw(4) << j            << "  "
         << setw(6) << g[i*ysize+j] << "\n";
  }

  delete [] g;

  return;
}
//****************************************************************************80

void test03 ( )  

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PGMA_WRITE.
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
{
# define NGRAY 11

  string output_name = "pgma_io_prb_03.ascii.pgm";
  int *g;
  double gray[NGRAY] = { 
    0.000, 0.291, 0.434, 0.540, 0.629,
    0.706, 0.774, 0.837, 0.895, 0.949,
    1.000 };
  int i;
  int j;
  int k;
  int xsize = 300;
  int ysize = 300;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  PGMA_WRITE writes an ASCII PGM file.\n";
  cout << "\n";
  cout << "  In this example, we make a sort of grayscale\n";
  cout << "  checkerboard.\n";

  g = new int[xsize*ysize];

  for ( i = 0; i < xsize; i++ )
  {
    for ( j = 0; j < ysize; j++ )
    {
      k = ( i + j ) * NGRAY / i4_min ( xsize, ysize );
      k = k % NGRAY;
      g[i*ysize+j] = ( int ) ( 255.0E+00 * gray[k] );
    }
  }

  cout << "  Writing the file \"" << output_name << "\".\n";

  pgma_write ( output_name, xsize, ysize, g );

  cout << "\n";
  cout <<  "  PGMA_WRITE was successful.\n";

  delete [] g;

  return;
# undef NGRAY
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
