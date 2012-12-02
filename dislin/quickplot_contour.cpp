# include <cstdlib>
# include <iostream>
# include <cmath>

# include "dislin.hpp"

using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUICKPLOT_CONTOUR demonstrates the DISLIN quickplot command QPLCON.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2012
//
//  Reference:
//
//    Helmut Michels,
//    The Data Plotting Software DISLIN - version 10.4,
//    Shaker Media GmbH, January 2010,
//    ISBN13: 978-3-86858-517-9.
//
{
  int i;
  int j;
  int levels;
  int m = 100;
  int n = 100;
  float x;
  float y;
  float *zmat;

  cout << "\n";
  cout << "QUICKPLOT_SURFACE:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the DISLIN 'quickplot' command QPLCON\n";
  cout << "  to make a contour plot of data stored as a matrix.\n";
//
//  Set up the data.
//
  zmat = new float[m*n];

  for ( i = 0; i < m; i++ )
  {
    x = 1.6 * ( float ) ( i ) / ( float ) ( m - 1 );
    for ( j = 0; j < n; j++ )
    {
      y = 1.6 * ( float ) ( j ) / ( float ) ( n - 1 );
      zmat[i+j*m] = pow ( x * x - 1.0, 2 ) + pow ( y * y - 1.0, 2 );
    }
  }
//
//  Specify the format of the output file.
//
  metafl ( "png" );
//
//  Specify that if a file already exists of the given name,
//  the new data should overwrite the old.
//
  filmod ( "delete" );
//
//  Specify the name of the output graphics file.
//
  setfil ( "quickplot_contour.png" );
//
//  Choose the page size and orientation.
//
  setpag ( "usal" );
//
//  For PNG output, reverse the default black background to white.
//
  scrmod ( "reverse" );
//
//  Open DISLIN.
//
  disini ( );
//
//  Label the axes and the plot.
//
  name ( "<-- X -->", "X" );
  name ( "<-- Y -->", "Y" );
  titlin ( "Quick plot by QPLCON", 2 );
//
//  Draw the curve.
//
  levels = 20;
  qplcon ( zmat, m, n, levels );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Free memory.
//
  delete [] zmat;
//
//  Terminate.
//
  cout << "\n";
  cout << "QUICKPLOT_SURFACE:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
