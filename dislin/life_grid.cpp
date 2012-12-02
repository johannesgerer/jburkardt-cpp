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
//    LIFE_GRID uses DISLIN to draw a grid for the game of Life.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 April 2011
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
# define N 100

  char ctit[60];
  int i;
  int j;
  int nr;
  int nx;
  int ny;
  int pat;
  float r;
  float x;
  float xvec[2];
  float y;
  float yvec[2];

  cout << "\n";
  cout << "LIFE_GRID:\n";
  cout << "  C++ version:\n";
  cout << "  Use DISLIN routines to plot a grid for Life..\n";
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
  setfil ( "life_grid.png" );
//
//  Choose the page size and orientation.
//  'USA' is 2160 plot units wide and 2790 plot units high.
//  'P' requests PROFILE rather than LANDSCAPE orientation.
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
//
//  Use a color table, which is required if we want to do color graphics.
//  For this color table, in particular,
//  1 = black,
//  2 = red,
//  3 = green,
//  4 = blue.
//
  setvlt ( "small" );
//
//  Define the X and Y sizes of the axis system in plot units.
//
  axslen ( 1500, 1500 );
//
//  Specify how the lower X, left Y, upper X and right Y axes are labeled.
//
  setgrf ( "line", "line", "line", "line" );
//
//  Set the axis origin 330 plot units to the right, and 2000 plot units DOWN.
//
  axspos ( 330, 2000 );
//
//  Relate the physical coordinates to the axes.
//
  graf ( 0.0, 100.0, 0.0, 0.5, 0.0, 100.0, 0.0, 0.5 );
//
//  Draw 21 horizontal lines.
//
  for ( j = 0; j <= 100; j = j + 5 )
  {
    y = ( float ) ( j );
    xvec[0] = 0.0;
    xvec[1] = 100.0;
    yvec[0] = y;
    yvec[1] = y;
    curve ( xvec, yvec, 2 );
  }
//
//  Draw 21 vertical lines.
//
  for ( i = 0; i <= 100; i = i + 5 )
  {
    x = ( float ) ( i );
    xvec[0] = x;
    xvec[1] = x;
    yvec[0] = 0.0;
    yvec[1] = 100.0;
    curve ( xvec, yvec, 2 );
  }
//
//  Select the shading pattern.
//
  pat = 16;
  shdpat ( pat );
//
//  Select color 3 (green) from the color table.
//
  setclr ( 3 );
//
//  Draw one circle near the origin.
//
  x = 2.5;
  y = 2.5;
  r = 2.0;
  rlcirc ( x, y, r );
//
//  Select color 2 (red).
//
  setclr ( 2 );
//
//  Draw a glider.
//
  x = 7.5;
  y = 37.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 12.5;
  y = 37.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 12.5;
  y = 47.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 17.5;
  y = 37.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 17.5;
  y = 42.5;
  r = 2.0;
  rlcirc ( x, y, r );
//
//  Select color 4 (blue)
//
  setclr ( 4 );
//
//  Select open shading pattern.
//
  pat = 0;
  shdpat ( pat );
//
//  Draw three open circles.
//
  x = 62.5;
  y = 62.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 67.5;
  y = 57.5;
  r = 2.0;
  rlcirc ( x, y, r );

  x = 72.5;
  y = 52.5;
  r = 2.0;
  rlcirc ( x, y, r );
//
//  Select character height in plot units.
//
  height ( 50 );
//
//  Select color 1 (black) from the color table.
//
  setclr ( 1 );
//
//  Define axis system titles.
//
  titlin ( "Grid for Game of Life", 1 );
//
//  Draw the title.
//
  title ( );
//
//  End this plot.
//
  endgrf ( );
//
//  Close DISLIN.
//
  disfin ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LIFE_GRID:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
