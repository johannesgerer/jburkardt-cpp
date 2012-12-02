# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

# include "image_edge.hpp"

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for IMAGE_EDGE_PRB.
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
  int *e;
  int *g;
  int *g_histo;
  int g_max;
  int i;
  string input_filename = "coins.ascii.pgm";
  ifstream input_unit;
  int m;
  int n;
  string output_filename = "coin_edges.ascii.pbm";

  timestamp ( );
  cout << "\n";
  cout << "IMAGE_EDGE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the IMAGE_EDGE library.\n";
  cout << "\n";
  cout << "  Demonstrate the NEWS stencil for edge detection\n";
  cout << "  in images.\n";

  cout << "\n";
  cout << "  The input file is \"" << input_filename << "\".\n";
//
//  Open the input file and read the data.
//
  input_unit.open ( input_filename.c_str ( ) );

  if ( !input_unit )
  {
    cerr << "\n";
    cerr << "IMAGE_EDGE_PRB - Fatal error!\n";
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

  g_histo = i4mat_histogram ( m, n, g, 255 );

  cout << "\n";
  cout << " Gray     Count\n";
  cout << "\n";
  for ( i = 0; i <= 255; i++ )
  {
    cout << "  " << setw(3) << i
         << "  " << setw(8) << g_histo[i] << "\n";
  }

  delete [] g_histo;

  e = news ( m, n, g );
//
//  Write the edge information as a portable BIT map (0/1).
//
  pbma_write ( output_filename, m, n, e );

  cout << "\n";
  cout << "  Wrote edge information to \"" << output_filename << "\".\n";

  delete [] e;
  delete [] g;
//
//  Terminate.
//
  cout << "\n";
  cout << "IMAGE_EDGE_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
