# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "image_denoise.hpp"

int main ( int argc, char *argv[] );
void test01 ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for IMAGE_DENOISE_PRB.
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
{
  timestamp ( );
  cout << "\n";
  cout << "IMAGE_DENOISE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the IMAGE_DENOISE library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "IMAGE_DENOISE_PRB\n";
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
//    TEST01 tests GRAY_MEDIAN_NEWS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  int *g;
  int g_max;
  int *g2;
  string input_filename = "glassware_noisy.ascii.pgm";
  ifstream input_unit;
  int m;
  int n;
  string output_filename = "glassware_median_news.ascii.pgm";

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  GRAY_MEDIAN_NEWS uses a NEWS median filter \n";
  cout << "  on a noisy grayscale image.\n";

  cout << "\n";
  cout << "  The input file is \"" << input_filename << "\".\n";
//
//  Open the input file and read the data.
//
  input_unit.open ( input_filename.c_str ( ) );

  if ( !input_unit )
  {
    cerr << "\n";
    cerr << "TEST01 - Fatal error!\n";
    cerr << "  Could not open the file \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  pgma_read_header ( input_unit, &m, &n, &g_max );

  cout << "\n";
  cout << "  Number of rows =          " << m << "\n";
  cout << "  Number of columns =       " << n << "\n";
  cout << "  Maximum pixel intensity = " << g_max << "\n";

  g = new int[m*n];

  pgma_read_data ( input_unit, m, n, g );

  input_unit.close ( );

  g2 = gray_median_news ( m, n, g );
//
//  Write the denoised images.
//
  pgma_write ( output_filename, m, n, g2 );

  cout << "\n";
  cout << "  Wrote denoised image to \"" << output_filename << "\".\n";

  delete [] g;
  delete [] g2;

  return;
}
