# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "dislin.hpp"

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN demonstrates the creation of a map plot.
//
//  Modified:
//
//    09 April 2011
//
//  Reference:
//
//    Helmut Michels,
//    The Data Plotting Software DISLIN - version 10.4,
//    Shaker Media GmbH, January 2010,
//    ISBN13: 978-3-86858-517-9.
//
{
  cout << "\n";
  cout << "DISLIN_EX13:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the creation of a map plot.\n";
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
  setfil ( "dislin_ex13.png" );
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
//  Plot a border around the page.
//
  pagera ( );
//
//  Use the COMPLEX font.
//
  complx ( );

  frame ( 3 );
  axspos ( 400, 1850 );
  axslen ( 2400, 1400 );

  name ( "Longitude", "x" );
  name ( "Latitude", "y" );
  titlin ( "World Coastlines and Lakes", 3 );

  labels ( "map", "xy" );
  grafmp ( -180.0, 180.0, -180.0, 90.0, -90.0, 90.0, -90.0, 30.0 );

  gridmp ( 1, 1 );
  color ( "green" );
  world ( );
  color ( "fore" );

  height ( 50 );
  title ( );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX13:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
