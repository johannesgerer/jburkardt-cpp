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
//    QUICKPLOT_PIE demonstrates the DISLIN quickplot command QPLPIE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 May 2012
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
  int n = 5;
  float xray[5] = { 10.0, 20.0, 15.0, 5.0, 50.0 };

  cout << "\n";
  cout << "QUICKPLOT_PIE:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the DISLIN 'quickplot' command QPLOT\n";
  cout << "  to plot a curve.\n";
  cout << "\n";
  cout << "  Here, we plot 10 percent luck, 20 percent skill,\n";
  cout << "  15 percent concentrated power of will, 5 percent pleasure,\n";
  cout << "  50 percent pain.\n";
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
  setfil ( "quickplot_pie.png" );
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
  titlin ( "Quick plot by QPLPIE", 2 );
//
//  Draw the curve.
//
  qplpie ( xray, n );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUICKPLOT_PIE:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
