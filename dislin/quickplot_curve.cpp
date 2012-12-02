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
//    QUICKPLOT_CURVE demonstrates the DISLIN quickplot command QPLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 April 2011
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
  int n = 100;
  float pi = 3.1415926;
  float *xray;
  float *yray;

  cout << "\n";
  cout << "QUICKPLOT_CURVE:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the DISLIN 'quickplot' command QPLOT\n";
  cout << "  to plot a curve.\n";
//
//  Set up the X and Y data for the plot.
//
  xray = new float[n];
  yray = new float[n];

  for ( i = 0; i < n; i++ )
  {
    xray[i] = ( float ) ( i ) * 360.0 / ( float ) ( n - 1 );
  }
  for ( i = 0; i < n; i++ )
  {
    yray[i] = sin ( pi * xray[i] / 180.0 );
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
  setfil ( "quickplot_curve.png" );
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
  name ( "<-- Angle in Degrees -->", "X" );
  name ( "<-- Sine (angle) -->", "Y" );
  titlin ( "Quick plot by QPLOT", 2 );
//
//  Draw the curve.
//
  qplot ( xray, yray, n );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Free memory.
//
  delete [] xray;
  delete [] yray;
//
//  Terminate.
//
  cout << "\n";
  cout << "QUICKPLOT_CURVE:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
