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
//    QUICKPLOT_BAR demonstrates the DISLIN quickplot command QPLBAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Helmut Michels,
//    The Data Plotting Software DISLIN - version 10.4,
//    Shaker Media GmbH, January 2010,
//    ISBN13: 978-3-86858-517-9.
//
{
# define N 14

  int i;
  int n = N;
  float xray[N] = {
     1.0, 15.0, 38.0, 22.0, 16.0, 
    16.0, 26.0, 55.0, 50.0, 40.0, 
    16.0,  3.0,  0.0,  1.0 };
    
  cout << "\n";
  cout << "QUICKPLOT_BAR:\n";
  cout << "  C++ version\n";
  cout << "  Demonstrate the DISLIN \"quickplot\" command QPLBAR\n";
  cout << "  to plot a bar chart.\n";
//
//  Specify the format of the output file.
//
  metafl ( "png" );
//
//  Indicate that new data overwrites old data.
//
  filmod ( "delete" );
//
//  Specify the name of the output graphics file.
//
  setfil ( "quickplot_bar.png" );
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
  name ( "<-- Minutes -->", "X" );
  name ( "<-- Frequency -->", "Y" );
  titlin ( "Quick plot by QPLBAR", 2 );
//
//  Draw the curve.
//
  qplbar ( xray, n );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUICKPLOT_BAR:\n";
  cout << "  Normal end of execution.\n";

  return 0;
# undef N
}
