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
//    MAIN demonstrates the use of the interpolation options when using CURVE.
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
  static char *cpol[6] = { "SPLINE", "STEM", "BARS", "STAIRS", "STEP", "LINEAR" };
  static char *ctit    = "Interpolation Methods";
  int i;
  int nx;
  int ny;
  int nya = 2700;
  static float x[] = {
     0.0,  1.0,  3.0,  4.5,  6.0,  8.0,  9.0, 11.0, 12.0, 12.5,
    13.0, 15.0, 16.0, 17.0, 19.0, 20.0 };
  static float y[] = {
     2.0,  4.0,  4.5,  3.0,  1.0,  7.0,  2.0,  3.0,  5.0,  2.0, 
     2.5,  2.0,  4.0,  6.0,  5.5,  4.0 };

  cout << "\n";
  cout << "DISLIN_EX04:\n";
  cout << "  C++ version:\n";
  cout << "  Demonstrate the use of the various interpolation options\n";
  cout << "  when using CURVE to plot (X,Y) data.\n";
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
  setfil ( "dislin_ex04.png" );
//
//  Choose the page size and orientation.
//
  setpag ( "usap" );
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
  incmrk ( 1 );
  hsymbl ( 25 );
  titlin ( ctit, 1 );
  axslen ( 1500, 350 );

  setgrf ( "line", "line", "line", "line" );

  for ( i = 0; i < 6; i++ )
  {
    axspos ( 350, nya-i*350 );
    polcrv ( cpol[i] );
    marker ( 0 );

    graf ( 0.0, 20.0, 0.0, 5.0, 0.0, 10.0, 0.0, 5.0 );
    nx = nxposn ( 1.0 );
    ny = nyposn ( 8.0 );
    messag ( cpol[i], nx, ny );
    curve ( x, y, 16 );

    if ( i == 5 )
    {
      height ( 50 );
      title ( );
    }
    endgrf ( );
  }
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DISLIN_EX04:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
